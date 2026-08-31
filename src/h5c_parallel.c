/*
 * Parallel (MPI-IO) implementation. Compiled only with H5C_ENABLE_PARALLEL.
 *
 * Layout and split axis are documented in include/h5c/h5c_mpi.h. The three
 * places where this deliberately differs from h5fortran are marked with
 * "DELIBERATE:" comments below.
 */
#include "h5c_internal.h"

#include "h5c/h5c_mpi.h"

#include <stdlib.h>
#include <string.h>

/* Name of the payload dataset inside the group written for a path. */
#define DATA_NAME "data"

/* ------------------------------------------------------------------ */
/* collective agreement                                                */
/* ------------------------------------------------------------------ */

/*
 * DELIBERATE: h5fortran lets a rank that fails validation return early while
 * the others walk into a collective HDF5 call, which deadlocks. Every status
 * here is agreed upon first, so all ranks leave together with the same value.
 * The enum is append-only and ordered with H5C_OK == 0, so MAX picks a real
 * failure over success.
 */
static h5c_status_t agree(MPI_Comm comm, h5c_status_t local)
{
    int mine, worst;

    mine  = (int)local;
    worst = mine;
    if (MPI_Allreduce(&mine, &worst, 1, MPI_INT, MPI_MAX, comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI, "MPI_Allreduce failed agreeing on status");
    }
    if (worst == (int)H5C_OK) {
        return H5C_OK;
    }
    if (local != H5C_OK) {
        return local; /* keep this rank's own, more specific message */
    }
    return h5c__fail((h5c_status_t)worst,
                     "another rank reported '%s'; failing collectively",
                     h5c_status_string((h5c_status_t)worst));
}

/* Validates that `file` was opened through the parallel entry points. */
static h5c_status_t pfile_check(h5c_file_t *file)
{
    if (file == NULL || file->fid < 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "file handle is NULL or closed");
    }
    if (!file->parallel) {
        return h5c__fail(H5C_ERR_STATE,
                         "file was not opened with h5c_popen()");
    }
    return H5C_OK;
}

/*
 * Local, non-collective argument checks shared by pwrite and pread.
 *
 * ONE DIMENSION POLICY: the handle, path, rank range and dims pointer are
 * validated by the serial checker, zero extents included (a rank may own no
 * rows). The only rule added here is the one that is genuinely parallel;
 * cross-rank agreement of rank and dims[1..] is collective and lives in
 * agree_shape().
 */
static h5c_status_t pcheck_args(h5c_file_t *file, const char *path,
                                int rank, const size_t *dims)
{
    h5c_status_t st;

    if ((st = pfile_check(file)) != H5C_OK) {
        return st;
    }
    if ((st = h5c__check_common(file, path, rank, dims)) != H5C_OK) {
        return st;
    }
    if (rank < 1) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "rank %d is invalid for '%s': "
                         "a scalar has no axis to split", rank, path);
    }
    return H5C_OK;
}

/*
 * DELIBERATE: the non-split dimensions must be identical on every rank.
 * Checked with MPI_Allreduce MIN/MAX (as h5fortran does) so the verdict is
 * already the same everywhere and no further agreement is needed.
 */
static h5c_status_t agree_shape(MPI_Comm comm, const char *path,
                                int rank, const size_t *dims)
{
    int64_t mine[H5C_MAX_RANK + 1];
    int64_t lo[H5C_MAX_RANK + 1];
    int64_t hi[H5C_MAX_RANK + 1];
    int     n, i;

    /*
     * A fixed count keeps the reduction well defined even when the ranks
     * disagree about `rank` itself. Element 0 is the rank, so that mismatch
     * is reported first; the unused tail is zero on every process.
     */
    n = H5C_MAX_RANK + 1;
    for (i = 0; i < n; i++) {
        mine[i] = 0;
    }
    mine[0] = (int64_t)rank;
    for (i = 1; i < rank; i++) {
        mine[i] = (int64_t)dims[i];
    }

    if (MPI_Allreduce(mine, lo, n, MPI_INT64_T, MPI_MIN, comm) != MPI_SUCCESS ||
        MPI_Allreduce(mine, hi, n, MPI_INT64_T, MPI_MAX, comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI,
                         "MPI_Allreduce failed checking the shape of '%s'",
                         path);
    }
    if (lo[0] != hi[0]) {
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         "'%s': rank differs across ranks (%ld..%ld)",
                         path, (long)lo[0], (long)hi[0]);
    }
    for (i = 1; i < n; i++) {
        if (lo[i] != hi[i]) {
            return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                             "'%s': dims[%d] differs across ranks "
                             "(%ld..%ld); only the split axis may differ",
                             path, i, (long)lo[i], (long)hi[i]);
        }
    }
    return H5C_OK;
}

/* Dataset transfer property list honouring the file's collective flag. */
static hid_t make_dxpl(const h5c_file_t *file)
{
    hid_t xfer;

    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (xfer < 0) {
        h5c__fail_hdf5((long)xfer, "H5Pcreate(H5P_DATASET_XFER) failed");
        return H5I_INVALID_HID;
    }
    if (H5Pset_dxpl_mpio(xfer, file->collective ? H5FD_MPIO_COLLECTIVE
                                                : H5FD_MPIO_INDEPENDENT) < 0) {
        H5Pclose(xfer);
        h5c__fail_hdf5(-1, "H5Pset_dxpl_mpio failed");
        return H5I_INVALID_HID;
    }
    return xfer;
}

/*
 * Selects this rank's block in `fsid` and builds a matching memory space.
 *
 * DELIBERATE: a local extent of 0 selects nothing at all (H5Sselect_none on
 * BOTH spaces) instead of a zero-length hyperslab, whose handling varies
 * between HDF5 and MPI-IO versions. The rank still joins the collective call.
 */
static h5c_status_t select_block(hid_t fsid, int rank, const size_t *dims,
                                 size_t offset, hid_t *msid_out)
{
    hsize_t start[H5C_MAX_RANK];
    hsize_t count[H5C_MAX_RANK];
    hsize_t mdims[H5C_MAX_RANK];
    hid_t   msid;
    int     i;

    *msid_out = H5I_INVALID_HID;

    for (i = 0; i < rank; i++) {
        start[i] = 0;
        count[i] = (hsize_t)dims[i];
        mdims[i] = (hsize_t)dims[i];
    }
    start[0] = (hsize_t)offset;
    mdims[0] = (dims[0] > 0) ? (hsize_t)dims[0] : 1;

    msid = H5Screate_simple(rank, mdims, NULL);
    if (msid < 0) {
        return h5c__fail_hdf5((long)msid, "cannot build the memory dataspace");
    }

    if (dims[0] == 0) {
        if (H5Sselect_none(fsid) < 0 || H5Sselect_none(msid) < 0) {
            H5Sclose(msid);
            return h5c__fail_hdf5(-1, "H5Sselect_none failed");
        }
    } else if (H5Sselect_hyperslab(fsid, H5S_SELECT_SET, start, NULL,
                                   count, NULL) < 0) {
        H5Sclose(msid);
        return h5c__fail_hdf5(-1, "H5Sselect_hyperslab failed");
    }

    *msid_out = msid;
    return H5C_OK;
}

/* ------------------------------------------------------------------ */
/* interleave tiling plan                                              */
/* ------------------------------------------------------------------ */

/*
 * Row-wise tiling for the interleaved entry points. The transfer loop below
 * is shared by pwrite_impl() and pread_impl(); the buffer arithmetic itself
 * comes from h5c_interleaved.c (h5c__tile_rows / h5c__pack_tile /
 * h5c__unpack_tile) so that serial and parallel cannot drift apart.
 *
 * `stage` holds `rows` whole rows and is NULL only when this rank owns none.
 * `ntiles` is AGREED ACROSS RANKS, see agree_tiles().
 */
typedef struct {
    void *const *comps;  /* borrowed; the write path const-casts into this */
    size_t       ncomp;
    size_t       esize;
    char        *stage;
    size_t       rows;   /* rows per tile; 0 when this rank owns no rows */
    long long    ntiles; /* identical on every rank */
} tile_plan_t;

/*
 * DELIBERATE: the tile count is agreed with MPI_Allreduce(MAX). Local row
 * counts differ between ranks, so a locally derived count would make ranks
 * issue different numbers of collective transfers and deadlock. Ranks that
 * run out of local rows before the agreed count is reached still enter every
 * call, with an empty selection (see select_tile). The floor of 1 keeps one
 * collective transfer even when no rank owns anything.
 */
static h5c_status_t agree_tiles(MPI_Comm comm, size_t n, size_t rows,
                                long long *ntiles)
{
    long long mine, most;

    mine = (rows > 0) ? (long long)((n + rows - 1) / rows) : 0;
    most = mine;
    if (MPI_Allreduce(&mine, &most, 1, MPI_LONG_LONG, MPI_MAX,
                      comm) != MPI_SUCCESS) {
        *ntiles = 1;
        return h5c__fail(H5C_ERR_MPI,
                         "MPI_Allreduce failed agreeing on the tile count");
    }
    *ntiles = (most < 1) ? 1 : most;
    return H5C_OK;
}

/*
 * Selects file rows [grow0, grow0 + rows) x ncomp and builds a matching
 * contiguous memory space. `rows == 0` selects nothing at all on both spaces,
 * exactly as select_block() does for an empty local block.
 */
static h5c_status_t select_tile(hid_t fsid, size_t grow0, size_t rows,
                                size_t ncomp, hid_t *msid_out)
{
    hsize_t start[2], count[2], mdims[2];
    hid_t   msid;

    *msid_out = H5I_INVALID_HID;

    mdims[0] = (rows > 0) ? (hsize_t)rows : 1;
    mdims[1] = (hsize_t)ncomp;
    msid = H5Screate_simple(2, mdims, NULL);
    if (msid < 0) {
        return h5c__fail_hdf5((long)msid, "cannot build the tile memory space");
    }

    if (rows == 0) {
        if (H5Sselect_none(fsid) < 0 || H5Sselect_none(msid) < 0) {
            H5Sclose(msid);
            return h5c__fail_hdf5(-1, "H5Sselect_none failed for a tile");
        }
    } else {
        start[0] = (hsize_t)grow0;
        start[1] = 0;
        count[0] = (hsize_t)rows;
        count[1] = (hsize_t)ncomp;
        if (H5Sselect_hyperslab(fsid, H5S_SELECT_SET, start, NULL,
                                count, NULL) < 0) {
            H5Sclose(msid);
            return h5c__fail_hdf5(-1, "cannot select tile rows");
        }
    }
    *msid_out = msid;
    return H5C_OK;
}

/*
 * Packs (writing) or unpacks (reading) `plan->ntiles` tiles, one collective
 * transfer each. A transfer error is remembered but does NOT leave the loop:
 * every rank must issue the same number of collective calls.
 */
static h5c_status_t transfer_tiles(hid_t did, hid_t fsid, hid_t mtype,
                                   hid_t xfer, const char *path, size_t n,
                                   size_t offset, const tile_plan_t *plan,
                                   int writing)
{
    h5c_status_t st = H5C_OK;
    char         dummy = 0;
    long long    t;

    for (t = 0; t < plan->ntiles; t++) {
        size_t row0 = (size_t)t * plan->rows;
        size_t rows = 0;
        char  *buf;
        hid_t  msid;

        if (plan->rows > 0 && row0 < n) {
            rows = (n - row0 < plan->rows) ? (n - row0) : plan->rows;
        }
        {
            h5c_status_t sel = select_tile(fsid, offset + row0, rows,
                                          plan->ncomp, &msid);
            if (sel != H5C_OK) {
                /* Unreachable in practice; nothing is left to select. */
                return sel;
            }
        }
        buf = (plan->stage != NULL) ? plan->stage : &dummy;

        if (writing) {
            h5c__pack_tile(buf, (const void *const *)plan->comps, plan->ncomp,
                           row0, rows, plan->esize);
            if (H5Dwrite(did, mtype, msid, fsid, xfer, buf) < 0 &&
                st == H5C_OK) {
                st = h5c__fail_hdf5(-1, "H5Dwrite failed for tile %lld of "
                                    "'%s/" DATA_NAME "'", t, path);
            }
        } else if (H5Dread(did, mtype, msid, fsid, xfer, buf) < 0) {
            if (st == H5C_OK) {
                st = h5c__fail_hdf5(-1, "H5Dread failed for tile %lld of "
                                    "'%s/" DATA_NAME "'", t, path);
            }
        } else {
            h5c__unpack_tile(buf, plan->comps, plan->ncomp, row0, rows,
                             plan->esize);
        }
        H5Sclose(msid);
    }
    return st;
}

/* ------------------------------------------------------------------ */
/* open / mode                                                         */
/* ------------------------------------------------------------------ */

h5c_status_t h5c_popen_comm(const char *path, h5c_mode_t mode,
                            MPI_Comm comm, MPI_Info info, h5c_file_t **out)
{
    h5c_status_t st;
    h5c_file_t  *file;
    hid_t        fapl, fid;

    if (out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_popen: out is NULL");
    }
    *out = NULL;
    if (path == NULL || path[0] == '\0') {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_popen: empty path");
    }
    if (comm == MPI_COMM_NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_popen: MPI_COMM_NULL");
    }
    if (mode != H5C_READ && mode != H5C_READWRITE && mode != H5C_TRUNCATE) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_popen: bad mode %d",
                         (int)mode);
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl < 0) {
        return h5c__fail_hdf5((long)fapl, "H5Pcreate(H5P_FILE_ACCESS) failed");
    }
    if (H5Pset_fapl_mpio(fapl, comm, info) < 0) {
        H5Pclose(fapl);
        return h5c__fail_hdf5(-1, "H5Pset_fapl_mpio failed for '%s'", path);
    }

    switch (mode) {
    case H5C_TRUNCATE:
        fid = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        break;
    case H5C_READ:
        fid = H5Fopen(path, H5F_ACC_RDONLY, fapl);
        break;
    default:
        fid = H5Fopen(path, H5F_ACC_RDWR, fapl);
        break;
    }
    H5Pclose(fapl);

    if (fid < 0) {
        return h5c__fail_hdf5((long)fid,
                              "cannot open '%s' in parallel (mode %d)",
                              path, (int)mode);
    }

    file = (h5c_file_t *)calloc(1, sizeof *file);
    if (file == NULL) {
        H5Fclose(fid);
        return h5c__fail(H5C_ERR_NOMEM, "h5c_popen: allocation failed");
    }
    file->fid        = fid;
    file->sticky     = H5C_OK;
    file->borrowed   = 0;
    file->readonly   = (mode == H5C_READ);
    file->parallel   = 1;
    file->collective = 1;  /* h5c never changes this implicitly */
    file->comm       = comm;  /* borrowed; valid until h5c_close() */

    *out = file;
    return H5C_OK;
}

h5c_status_t h5c_popen(const char *path, h5c_mode_t mode, h5c_file_t **out)
{
    return h5c_popen_comm(path, mode, MPI_COMM_WORLD, MPI_INFO_NULL, out);
}

h5c_status_t h5c_pset_collective(h5c_file_t *file, int collective)
{
    h5c_status_t st;

    if ((st = pfile_check(file)) != H5C_OK) {
        return h5c__record(file, st);
    }
    file->collective = collective ? 1 : 0;
    return H5C_OK;
}

int h5c_pis_collective(const h5c_file_t *file)
{
    if (file == NULL || !file->parallel) {
        return 0;
    }
    return file->collective;
}

MPI_Comm h5c_pcomm(const h5c_file_t *file)
{
    if (file == NULL || !file->parallel) {
        return MPI_COMM_NULL;
    }
    return file->comm;
}

/* ------------------------------------------------------------------ */
/* partition dataset                                                   */
/* ------------------------------------------------------------------ */

/*
 * Writes __partition__ collectively, matching h5fortran's split: rank 0
 * contributes entries 0 and 1, rank r contributes entry r+1. Every rank holds
 * the whole vector already, but writing it this way keeps every rank inside
 * the collective call.
 */
static h5c_status_t write_partition(hid_t gid, const int64_t *part,
                                    int me, int nprocs, hid_t xfer)
{
    hsize_t fdims[1], mdims[1], start[1], count[1];
    hid_t   fsid, msid, did;
    herr_t  err;

    fdims[0] = (hsize_t)(nprocs + 1);
    if (me == 0) {
        start[0] = 0;
        count[0] = 2;
    } else {
        start[0] = (hsize_t)(me + 1);
        count[0] = 1;
    }
    mdims[0] = count[0];

    fsid = H5Screate_simple(1, fdims, NULL);
    if (fsid < 0) {
        return h5c__fail_hdf5((long)fsid, "cannot build the partition space");
    }
    did = H5Dcreate2(gid, H5C_PARTITION_NAME, H5T_STD_I64LE, fsid,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (did < 0) {
        H5Sclose(fsid);
        return h5c__fail_hdf5((long)did,
                              "cannot create " H5C_PARTITION_NAME);
    }
    msid = H5Screate_simple(1, mdims, NULL);
    if (msid < 0) {
        H5Dclose(did);
        H5Sclose(fsid);
        return h5c__fail_hdf5((long)msid,
                              "cannot build the partition memory space");
    }
    if (H5Sselect_hyperslab(fsid, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) {
        H5Sclose(msid);
        H5Dclose(did);
        H5Sclose(fsid);
        return h5c__fail_hdf5(-1, "cannot select the partition slab");
    }

    err = H5Dwrite(did, H5T_NATIVE_INT64, msid, fsid, xfer,
                   (me == 0) ? &part[0] : &part[me + 1]);

    H5Sclose(msid);
    H5Dclose(did);
    H5Sclose(fsid);

    if (err < 0) {
        return h5c__fail_hdf5(-1, "H5Dwrite failed for " H5C_PARTITION_NAME);
    }
    return H5C_OK;
}

/* Reads and validates __partition__. `data_rows` is the extent of data[0]. */
static h5c_status_t read_partition(hid_t gid, int nprocs, int64_t *part,
                                   size_t data_rows)
{
    h5c_dataset_info_t info;
    h5c_status_t       st;
    hid_t              did;
    int                r;

    did = H5Dopen2(gid, H5C_PARTITION_NAME, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND,
                         "no " H5C_PARTITION_NAME " beside the data");
    }
    if ((st = h5c__info_from_dset(did, &info)) != H5C_OK) {
        H5Dclose(did);
        return st;
    }
    if (info.rank != 1 || info.dims[0] != (size_t)(nprocs + 1)) {
        H5Dclose(did);
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         H5C_PARTITION_NAME " has %lu entries, expected %d "
                         "(one per rank plus one)",
                         (unsigned long)((info.rank == 1) ? info.dims[0] : 0),
                         nprocs + 1);
    }
    /* Independent metadata-sized read; every rank needs the whole vector. */
    if (H5Dread(did, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, part) < 0) {
        H5Dclose(did);
        return h5c__fail_hdf5(-1, "H5Dread failed for " H5C_PARTITION_NAME);
    }
    H5Dclose(did);

    if (part[0] != 0) {
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         H5C_PARTITION_NAME "[0] is %ld, expected 0",
                         (long)part[0]);
    }
    for (r = 0; r < nprocs; r++) {
        if (part[r + 1] < part[r]) {
            return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                             H5C_PARTITION_NAME " decreases at %d "
                             "(%ld then %ld)", r + 1,
                             (long)part[r], (long)part[r + 1]);
        }
    }
    if ((size_t)part[nprocs] != data_rows) {
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         H5C_PARTITION_NAME " ends at %ld but data has %lu "
                         "rows along the split axis",
                         (long)part[nprocs], (unsigned long)data_rows);
    }
    return H5C_OK;
}

/* ------------------------------------------------------------------ */
/* write                                                               */
/* ------------------------------------------------------------------ */

/*
 * Writes this rank's block. With `pack == NULL` the caller's contiguous `buf`
 * is written in one collective transfer; with a plan, the data is gathered
 * from the plan's components tile by tile (see transfer_tiles). Everything
 * else - validation, agreement, group and __partition__ handling - is shared.
 */
static h5c_status_t pwrite_impl(h5c_file_t *file, const char *path,
                                const void *buf, h5c_type_t type,
                                int rank, const size_t *dims, unsigned flags,
                                const tile_plan_t *pack)
{
    h5c_status_t st;
    MPI_Comm     comm;
    int64_t     *part = NULL;
    int64_t      nlocal;
    hsize_t      fdims[H5C_MAX_RANK];
    hid_t        ftype, mtype;
    hid_t        gid = H5I_INVALID_HID, did = H5I_INVALID_HID;
    hid_t        fsid = H5I_INVALID_HID, msid = H5I_INVALID_HID;
    hid_t        xfer = H5I_INVALID_HID;
    size_t       offset;
    int          me, nprocs, r, i, existed;
    char         dummy = 0;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }
    if (file == NULL || !file->parallel) {
        /* No usable communicator, so nothing can be agreed. */
        return pfile_check(file);
    }
    comm = file->comm;

    /* --- local validation, agreed before any HDF5 call ------------- */
    st = pcheck_args(file, path, rank, dims);
    if (st == H5C_OK && pack == NULL && buf == NULL && dims[0] > 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "buffer is NULL for '%s'", path);
    }
    if (st == H5C_OK && file->readonly) {
        st = h5c__fail(H5C_ERR_STATE,
                       "file is open read-only, cannot write '%s'", path);
    }
    if (st == H5C_OK) {
        ftype = h5c__file_type(type);
        mtype = h5c__mem_type(type);
        if (ftype == H5I_INVALID_HID || mtype == H5I_INVALID_HID) {
            st = h5c__fail(H5C_ERR_INVALID_ARG,
                           "type %d has no numeric mapping "
                           "(parallel string I/O is not supported)",
                           (int)type);
        }
    } else {
        ftype = mtype = H5I_INVALID_HID;
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        return st;
    }
    if ((st = agree_shape(comm, path, rank, dims)) != H5C_OK) {
        return st;  /* identical on every rank already */
    }

    if (MPI_Comm_rank(comm, &me) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &nprocs) != MPI_SUCCESS) {
        return agree(comm, h5c__fail(H5C_ERR_MPI, "MPI_Comm_rank/size failed"));
    }

    /* --- partition from the local extents -------------------------- */
    part = (int64_t *)calloc((size_t)nprocs + 1, sizeof *part);
    if (part == NULL) {
        return agree(comm, h5c__fail(H5C_ERR_NOMEM,
                                     "cannot allocate the partition vector"));
    }
    nlocal = (int64_t)dims[0];
    if (MPI_Allgather(&nlocal, 1, MPI_INT64_T, part + 1, 1, MPI_INT64_T,
                      comm) != MPI_SUCCESS) {
        free(part);
        return agree(comm, h5c__fail(H5C_ERR_MPI, "MPI_Allgather failed"));
    }
    part[0] = 0;
    for (r = 0; r < nprocs; r++) {
        part[r + 1] += part[r];
    }
    offset = (size_t)part[me];

    /* --- group, replacing an existing one only when asked ---------- */
    existed = h5c_exists(file, path);
    if (existed && (flags & H5C_WRITE_REPLACE)) {
        if (H5Ldelete(file->fid, path, H5P_DEFAULT) < 0) {
            st = h5c__fail_hdf5(-1, "cannot replace existing '%s'", path);
        }
        existed = 0;
    } else if (existed) {
        st = h5c__fail(H5C_ERR_EXISTS,
                       "'%s' already exists; pass H5C_WRITE_REPLACE "
                       "to overwrite it", path);
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        free(part);
        return st;
    }

    gid = H5Gcreate2(file->fid, path, h5c__lcpl(), H5P_DEFAULT, H5P_DEFAULT);
    if (gid < 0) {
        free(part);
        return h5c__fail_hdf5((long)gid, "cannot create the group '%s'", path);
    }

    xfer = make_dxpl(file);
    if (xfer == H5I_INVALID_HID) {
        H5Gclose(gid);
        free(part);
        return H5C_ERR_HDF5;
    }

    /* --- data ------------------------------------------------------ */
    fdims[0] = (hsize_t)part[nprocs];
    for (i = 1; i < rank; i++) {
        fdims[i] = (hsize_t)dims[i];
    }
    fsid = H5Screate_simple(rank, fdims, NULL);
    if (fsid < 0) {
        st = h5c__fail_hdf5((long)fsid, "cannot build the file dataspace");
        goto done;
    }
    did = H5Dcreate2(gid, DATA_NAME, ftype, fsid,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (did < 0) {
        st = h5c__fail_hdf5((long)did, "cannot create '%s/" DATA_NAME "'",
                            path);
        goto done;
    }
    if (pack != NULL) {
        st = transfer_tiles(did, fsid, mtype, xfer, path, dims[0], offset,
                            pack, 1);
    } else if ((st = select_block(fsid, rank, dims, offset, &msid)) == H5C_OK) {
        if (H5Dwrite(did, mtype, msid, fsid, xfer,
                     (buf != NULL) ? buf : (const void *)&dummy) < 0) {
            st = h5c__fail_hdf5(-1, "H5Dwrite failed for '%s/" DATA_NAME "'",
                                path);
        }
    }
    /* __partition__ is collective too, so agree before entering it. */
    if ((st = agree(comm, st)) != H5C_OK) {
        goto done;
    }

    st = write_partition(gid, part, me, nprocs, xfer);

done:
    if (msid >= 0) { H5Sclose(msid); }
    if (did  >= 0) { H5Dclose(did);  }
    if (fsid >= 0) { H5Sclose(fsid); }
    if (xfer >= 0) { H5Pclose(xfer); }
    if (gid  >= 0) { H5Gclose(gid);  }
    free(part);
    return st;
}

h5c_status_t h5c_pwrite(h5c_file_t *file, const char *path, const void *buf,
                        h5c_type_t type, int rank, const size_t *dims,
                        unsigned flags)
{
    return h5c__record(file, pwrite_impl(file, path, buf, type, rank, dims,
                                        flags, NULL));
}

/* ------------------------------------------------------------------ */
/* read                                                                */
/* ------------------------------------------------------------------ */

/*
 * Opens path/data, reads and validates __partition__, and reports this rank's
 * block. On success *gid_out and *did_out are open and owned by the caller.
 * The returned status is already agreed across the communicator.
 */
static h5c_status_t open_for_read(h5c_file_t *file, const char *path,
                                  int me, int nprocs,
                                  hid_t *gid_out, hid_t *did_out,
                                  h5c_dataset_info_t *info,
                                  int64_t *part, size_t *offset,
                                  size_t *nlocal)
{
    h5c_status_t st = H5C_OK;
    hid_t        gid = H5I_INVALID_HID, did = H5I_INVALID_HID;

    *gid_out = H5I_INVALID_HID;
    *did_out = H5I_INVALID_HID;

    gid = H5Gopen2(file->fid, path, H5P_DEFAULT);
    if (gid < 0) {
        st = h5c__fail(H5C_ERR_NOT_FOUND,
                       "no distributed dataset group at '%s'", path);
    }
    if (st == H5C_OK) {
        did = H5Dopen2(gid, DATA_NAME, H5P_DEFAULT);
        if (did < 0) {
            st = h5c__fail(H5C_ERR_NOT_FOUND, "no '%s/" DATA_NAME "'", path);
        }
    }
    if (st == H5C_OK) {
        st = h5c__info_from_dset(did, info);
    }
    if (st == H5C_OK && info->rank < 1) {
        st = h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                       "'%s/" DATA_NAME "' is a scalar; nothing to split",
                       path);
    }
    if (st == H5C_OK) {
        st = read_partition(gid, nprocs, part, info->dims[0]);
    }
    if (st == H5C_OK) {
        *offset = (size_t)part[me];
        *nlocal = (size_t)(part[me + 1] - part[me]);
    }

    /* Agree BEFORE the caller enters the collective transfer. */
    st = agree(file->comm, st);
    if (st != H5C_OK) {
        if (did >= 0) { H5Dclose(did); }
        if (gid >= 0) { H5Gclose(gid); }
        return st;
    }
    *gid_out = gid;
    *did_out = did;
    return H5C_OK;
}

/* Reads this rank's block; `unpack` mirrors pwrite_impl's `pack`. */
static h5c_status_t pread_impl(h5c_file_t *file, const char *path, void *buf,
                               h5c_type_t type, int rank, const size_t *dims,
                               const tile_plan_t *unpack)
{
    h5c_status_t       st;
    h5c_dataset_info_t info;
    MPI_Comm           comm;
    int64_t           *part = NULL;
    hid_t              mtype;
    hid_t              gid = H5I_INVALID_HID, did = H5I_INVALID_HID;
    hid_t              fsid = H5I_INVALID_HID, msid = H5I_INVALID_HID;
    hid_t              xfer = H5I_INVALID_HID;
    size_t             offset = 0, nlocal = 0;
    int                me, nprocs, i;
    char               dummy = 0;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }
    if (file == NULL || !file->parallel) {
        return pfile_check(file);
    }
    comm = file->comm;

    st = pcheck_args(file, path, rank, dims);
    if (st == H5C_OK && unpack == NULL && buf == NULL && dims[0] > 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "buffer is NULL for '%s'", path);
    }
    if (st == H5C_OK) {
        mtype = h5c__mem_type(type);
        if (mtype == H5I_INVALID_HID) {
            st = h5c__fail(H5C_ERR_INVALID_ARG,
                           "type %d has no numeric mapping "
                           "(parallel string I/O is not supported)",
                           (int)type);
        }
    } else {
        mtype = H5I_INVALID_HID;
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        return st;
    }
    if ((st = agree_shape(comm, path, rank, dims)) != H5C_OK) {
        return st;
    }

    if (MPI_Comm_rank(comm, &me) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &nprocs) != MPI_SUCCESS) {
        return agree(comm, h5c__fail(H5C_ERR_MPI, "MPI_Comm_rank/size failed"));
    }
    part = (int64_t *)calloc((size_t)nprocs + 1, sizeof *part);
    if (part == NULL) {
        return agree(comm, h5c__fail(H5C_ERR_NOMEM,
                                     "cannot allocate the partition vector"));
    }

    st = open_for_read(file, path, me, nprocs, &gid, &did,
                       &info, part, &offset, &nlocal);
    if (st != H5C_OK) {
        free(part);
        return st;
    }

    /* The caller's block must be exactly what __partition__ assigns it. */
    if (info.rank != rank) {
        st = h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                       "'%s/" DATA_NAME "' has rank %d, expected %d",
                       path, info.rank, rank);
    } else if (dims[0] != nlocal) {
        st = h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                       "'%s': rank %d owns %lu rows per " H5C_PARTITION_NAME
                       ", but dims[0] is %lu",
                       path, me, (unsigned long)nlocal,
                       (unsigned long)dims[0]);
    } else {
        for (i = 1; i < rank; i++) {
            if (info.dims[i] != dims[i]) {
                st = h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                               "'%s/" DATA_NAME "' dims[%d] is %lu, "
                               "expected %lu", path, i,
                               (unsigned long)info.dims[i],
                               (unsigned long)dims[i]);
                break;
            }
        }
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        goto done;
    }

    xfer = make_dxpl(file);
    if (xfer == H5I_INVALID_HID) {
        st = H5C_ERR_HDF5;
        goto done;
    }
    fsid = H5Dget_space(did);
    if (fsid < 0) {
        st = h5c__fail_hdf5((long)fsid, "H5Dget_space failed for '%s'", path);
        goto done;
    }
    if (unpack != NULL) {
        st = transfer_tiles(did, fsid, mtype, xfer, path, dims[0], offset,
                            unpack, 0);
    } else if ((st = select_block(fsid, rank, dims, offset, &msid)) == H5C_OK) {
        if (H5Dread(did, mtype, msid, fsid, xfer,
                    (buf != NULL) ? buf : (void *)&dummy) < 0) {
            st = h5c__fail_hdf5(-1, "H5Dread failed for '%s/" DATA_NAME "'",
                                path);
        }
    }

done:
    if (msid >= 0) { H5Sclose(msid); }
    if (fsid >= 0) { H5Sclose(fsid); }
    if (xfer >= 0) { H5Pclose(xfer); }
    if (did  >= 0) { H5Dclose(did);  }
    if (gid  >= 0) { H5Gclose(gid);  }
    free(part);
    return st;
}

h5c_status_t h5c_pread(h5c_file_t *file, const char *path, void *buf,
                       h5c_type_t type, int rank, const size_t *dims)
{
    return h5c__record(file, pread_impl(file, path, buf, type, rank, dims,
                                       NULL));
}

/* ------------------------------------------------------------------ */
/* shape query                                                         */
/* ------------------------------------------------------------------ */

static h5c_status_t pinfo_impl(h5c_file_t *file, const char *path,
                               h5c_dataset_info_t *local,
                               h5c_dataset_info_t *global)
{
    h5c_status_t       st;
    h5c_dataset_info_t info;
    MPI_Comm           comm;
    int64_t           *part = NULL;
    hid_t              gid, did;
    size_t             offset = 0, nlocal = 0;
    int                me, nprocs, i;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }
    if (file == NULL || !file->parallel) {
        return pfile_check(file);
    }
    comm = file->comm;

    st = h5c__check_common(file, path, 0, NULL);
    if (st == H5C_OK && local == NULL && global == NULL) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "h5c_pdataset_info: both outputs are NULL");
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        return st;
    }

    if (MPI_Comm_rank(comm, &me) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &nprocs) != MPI_SUCCESS) {
        return agree(comm, h5c__fail(H5C_ERR_MPI, "MPI_Comm_rank/size failed"));
    }
    part = (int64_t *)calloc((size_t)nprocs + 1, sizeof *part);
    if (part == NULL) {
        return agree(comm, h5c__fail(H5C_ERR_NOMEM,
                                     "cannot allocate the partition vector"));
    }

    st = open_for_read(file, path, me, nprocs, &gid, &did,
                       &info, part, &offset, &nlocal);
    free(part);
    if (st != H5C_OK) {
        return st;
    }
    H5Dclose(did);
    H5Gclose(gid);

    if (global != NULL) {
        *global = info;
    }
    if (local != NULL) {
        *local = info;
        local->dims[0] = nlocal;
        local->count   = nlocal;
        for (i = 1; i < info.rank; i++) {
            local->count *= info.dims[i];
        }
    }
    return H5C_OK;
}

h5c_status_t h5c_pdataset_info(h5c_file_t *file, const char *path,
                               h5c_dataset_info_t *local,
                               h5c_dataset_info_t *global)
{
    return h5c__record(file, pinfo_impl(file, path, local, global));
}

/* ------------------------------------------------------------------ */
/* layout accessors                                                    */
/* ------------------------------------------------------------------ */

/*
 * Reads and validates __partition__, then hands the caller the whole boundary
 * vector plus this rank's slice. Shared by the two public accessors so that
 * neither has to know how the layout is stored.
 *
 * On success `*part_out` is a malloc'd vector of `*nprocs_out + 1` entries and
 * belongs to the caller.
 */
static h5c_status_t players_impl(h5c_file_t *file, const char *path,
                                 h5c_status_t extra,
                                 int64_t **part_out, int *nprocs_out,
                                 size_t *offset_out, size_t *nlocal_out)
{
    h5c_status_t       st;
    h5c_dataset_info_t info;
    MPI_Comm           comm;
    int64_t           *part = NULL;
    hid_t              gid, did;
    int                me, nprocs;

    *part_out   = NULL;
    *nprocs_out = 0;
    *offset_out = 0;
    *nlocal_out = 0;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }
    if (file == NULL || !file->parallel) {
        return pfile_check(file);
    }
    comm = file->comm;

    st = h5c__check_common(file, path, 0, NULL);
    if (st == H5C_OK) {
        /*
         * The caller's own argument check rides along on this agreement. It
         * must NOT short-circuit earlier: a rank that returned before the
         * MPI_Allreduce would leave its peers waiting inside it.
         */
        st = extra;
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        return st;
    }

    if (MPI_Comm_rank(comm, &me) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &nprocs) != MPI_SUCCESS) {
        return agree(comm, h5c__fail(H5C_ERR_MPI, "MPI_Comm_rank/size failed"));
    }
    part = (int64_t *)calloc((size_t)nprocs + 1, sizeof *part);
    if (part == NULL) {
        return agree(comm, h5c__fail(H5C_ERR_NOMEM,
                                     "cannot allocate the partition vector"));
    }

    /* open_for_read validates __partition__ and agrees the status for us. */
    st = open_for_read(file, path, me, nprocs, &gid, &did,
                       &info, part, offset_out, nlocal_out);
    if (st != H5C_OK) {
        free(part);
        return st;
    }
    H5Dclose(did);
    H5Gclose(gid);

    *part_out   = part;
    *nprocs_out = nprocs;
    return H5C_OK;
}

static h5c_status_t poffset_impl(h5c_file_t *file, const char *path,
                                 size_t *offset, size_t *nlocal)
{
    h5c_status_t st;
    int64_t     *part = NULL;
    size_t       off = 0, n = 0;
    int          nprocs = 0;

    st = players_impl(file, path, H5C_OK, &part, &nprocs, &off, &n);
    if (st != H5C_OK) {
        return st;
    }
    free(part);

    if (offset != NULL) {
        *offset = off;
    }
    if (nlocal != NULL) {
        *nlocal = n;
    }
    return H5C_OK;
}

h5c_status_t h5c_poffset(h5c_file_t *file, const char *path,
                         size_t *offset, size_t *nlocal)
{
    return h5c__record(file, poffset_impl(file, path, offset, nlocal));
}

static h5c_status_t ppartition_impl(h5c_file_t *file, const char *path,
                                    int64_t *bounds, size_t capacity,
                                    size_t *count)
{
    h5c_status_t st;
    int64_t     *part = NULL;
    size_t       off = 0, n = 0, len;
    int          nprocs = 0, r;

    st = (bounds == NULL && count == NULL)
             ? h5c__fail(H5C_ERR_INVALID_ARG,
                         "h5c_ppartition: both outputs are NULL")
             : H5C_OK;

    st = players_impl(file, path, st, &part, &nprocs, &off, &n);
    if (st != H5C_OK) {
        return st;
    }
    len = (size_t)nprocs + 1;

    if (bounds != NULL && capacity < len) {
        st = h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                       "h5c_ppartition: capacity %lu is short of the %lu "
                       "boundaries this decomposition needs",
                       (unsigned long)capacity, (unsigned long)len);
    }
    if (st == H5C_OK && bounds != NULL) {
        for (r = 0; r < nprocs + 1; r++) {
            bounds[r] = part[r];
        }
    }
    free(part);

    if (st == H5C_OK && count != NULL) {
        *count = len;
    }
    /*
     * Keep the ranks in step even on this purely local check: a caller whose
     * capacity is wrong on one rank only should still see one shared verdict.
     * No collective HDF5 call follows, so this is the last chance to agree.
     */
    return agree(file->comm, st);
}

h5c_status_t h5c_ppartition(h5c_file_t *file, const char *path,
                            int64_t *bounds, size_t capacity, size_t *count)
{
    return h5c__record(file,
                       ppartition_impl(file, path, bounds, capacity, count));
}

/* ------------------------------------------------------------------ */
/* interleaved multi-component fields                                  */
/* ------------------------------------------------------------------ */

/* Local checks shared by the two interleaved entry points. */
static h5c_status_t pcheck_comps(const void *const *comps, size_t ncomp,
                                 size_t n, h5c_type_t type, size_t *esize)
{
    size_t c;

    *esize = h5c_type_size(type);
    if (*esize == 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d cannot be interleaved", (int)type);
    }
    if (ncomp == 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "ncomp is 0");
    }
    if (comps == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "comps is NULL");
    }
    if (n > 0) {
        for (c = 0; c < ncomp; c++) {
            if (comps[c] == NULL) {
                return h5c__fail(H5C_ERR_INVALID_ARG,
                                 "comps[%lu] is NULL", (unsigned long)c);
            }
        }
    }
    return H5C_OK;
}

h5c_status_t h5c_pwrite_interleaved(h5c_file_t *file, const char *path,
                                    const void *const *comps, size_t ncomp,
                                    size_t n, h5c_type_t type, unsigned flags)
{
    h5c_status_t st;
    MPI_Comm     comm;
    tile_plan_t  plan;
    size_t       dims[2], row_bytes;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return h5c__record(file, st);
    }
    if (file == NULL || !file->parallel) {
        return h5c__record(file, pfile_check(file));
    }
    comm = file->comm;

    memset(&plan, 0, sizeof plan);
    st = pcheck_comps(comps, ncomp, n, type, &plan.esize);
    if (st == H5C_OK) {
        plan.comps = (void *const *)comps;
        plan.ncomp = ncomp;
        row_bytes  = ncomp * plan.esize;
        plan.rows  = h5c__tile_rows(n, row_bytes);
        if (plan.rows > 0) {
            plan.stage = (char *)malloc(plan.rows * row_bytes);
            if (plan.stage == NULL) {
                st = h5c__fail(H5C_ERR_NOMEM,
                               "cannot allocate %lu bytes to pack '%s'",
                               (unsigned long)(plan.rows * row_bytes), path);
            }
        }
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        free(plan.stage);
        return h5c__record(file, st);
    }
    /* The tile count must be identical everywhere, see agree_tiles(). */
    if ((st = agree_tiles(comm, n, plan.rows, &plan.ntiles)) != H5C_OK) {
        free(plan.stage);
        return h5c__record(file, st);
    }

    /*
     * Pack, never stride: a strided collective write makes MPI-IO fill file
     * blocks partially and turns the write into read-modify-write. The
     * gathering itself happens tile by tile inside pwrite_impl().
     */
    dims[0] = n;
    dims[1] = ncomp;
    st = h5c__record(file,
                     pwrite_impl(file, path, NULL, type, 2, dims, flags,
                                 &plan));
    free(plan.stage);
    return st;
}

h5c_status_t h5c_pread_interleaved(h5c_file_t *file, const char *path,
                                   void *const *comps, size_t ncomp,
                                   size_t n, h5c_type_t type)
{
    h5c_status_t st;
    MPI_Comm     comm;
    tile_plan_t  plan;
    size_t       dims[2], row_bytes;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return h5c__record(file, st);
    }
    if (file == NULL || !file->parallel) {
        return h5c__record(file, pfile_check(file));
    }
    comm = file->comm;

    memset(&plan, 0, sizeof plan);
    st = pcheck_comps((const void *const *)comps, ncomp, n, type, &plan.esize);
    if (st == H5C_OK) {
        plan.comps = comps;
        plan.ncomp = ncomp;
        row_bytes  = ncomp * plan.esize;
        plan.rows  = h5c__tile_rows(n, row_bytes);
        if (plan.rows > 0) {
            plan.stage = (char *)malloc(plan.rows * row_bytes);
            if (plan.stage == NULL) {
                st = h5c__fail(H5C_ERR_NOMEM,
                               "cannot allocate %lu bytes to unpack '%s'",
                               (unsigned long)(plan.rows * row_bytes), path);
            }
        }
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        free(plan.stage);
        return h5c__record(file, st);
    }
    if ((st = agree_tiles(comm, n, plan.rows, &plan.ntiles)) != H5C_OK) {
        free(plan.stage);
        return h5c__record(file, st);
    }

    /*
     * All components are wanted, so one contiguous block read per tile beats
     * ncomp strided reads. h5c_read_component() is the strided path.
     */
    dims[0] = n;
    dims[1] = ncomp;
    st = h5c__record(file, pread_impl(file, path, NULL, type, 2, dims, &plan));
    free(plan.stage);
    return st;
}
