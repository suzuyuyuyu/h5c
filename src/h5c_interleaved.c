/*
 * Interleaved multi-component fields: [n, ncomp] datasets built from, or
 * scattered back into, `ncomp` separate component arrays of `n` elements.
 *
 * WRITING PACKS, READING STRIDES (see docs/h5c-h5cpp-design.md):
 * a strided collective write makes MPI-IO fill file blocks only partially and
 * forces read-modify-write, so writes gather into a contiguous staging buffer
 * and issue contiguous hyperslab writes. Reads carry no such penalty, so
 * h5c_read_component() selects a strided hyperslab and moves 1/ncomp of the
 * bytes.
 *
 * When n * ncomp * elemsize exceeds the pack limit the transfer is split into
 * ROW-WISE tiles: each tile still covers whole rows, so every file selection
 * stays contiguous while the staging buffer stays bounded.
 */
#include "h5c_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* pack limit                                                          */
/* ------------------------------------------------------------------ */

/*
 * Process-wide, deliberately a plain static: it is a tuning knob set once at
 * start-up, not per-call state. Concurrent set/get from several threads is not
 * synchronised.
 *
 * h5c_set_pack_limit(0) is CLAMPED (not rejected) to H5C_PACK_LIMIT_MIN, so a
 * caller passing a nonsense value still gets a working, merely tiny, tile.
 */
#define H5C_PACK_LIMIT_DEFAULT ((size_t)256u * 1024u * 1024u)
#define H5C_PACK_LIMIT_MIN     ((size_t)4096u)

static size_t g_pack_limit = H5C_PACK_LIMIT_DEFAULT;

void h5c_set_pack_limit(size_t bytes)
{
    g_pack_limit = (bytes < H5C_PACK_LIMIT_MIN) ? H5C_PACK_LIMIT_MIN : bytes;
}

size_t h5c_pack_limit(void)
{
    return g_pack_limit;
}

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

/*
 * Buffer-level helpers, declared in h5c_internal.h and shared with the
 * parallel entry points so that only one implementation of the interleave
 * arithmetic exists. The two public layers stay separate.
 */

/* Rows per tile: as many whole rows as the pack limit allows, at least 1. */
size_t h5c__tile_rows(size_t n, size_t row_bytes)
{
    size_t rows = g_pack_limit / row_bytes;

    if (rows == 0) {
        rows = 1;
    }
    return (rows > n) ? n : rows;
}

/* Gathers rows [row0, row0+rows) of every component into `dst`. */
void h5c__pack_tile(char *dst, const void *const *comps, size_t ncomp,
                    size_t row0, size_t rows, size_t esize)
{
    size_t c, r;

    for (c = 0; c < ncomp; c++) {
        const char *src = (const char *)comps[c] + row0 * esize;
        char       *out = dst + c * esize;
        for (r = 0; r < rows; r++) {
            memcpy(out, src, esize);
            src += esize;
            out += ncomp * esize;
        }
    }
}

/* Scatters rows [row0, row0+rows) of `src` back into the components. */
void h5c__unpack_tile(const char *src, void *const *comps, size_t ncomp,
                      size_t row0, size_t rows, size_t esize)
{
    size_t c, r;

    for (c = 0; c < ncomp; c++) {
        const char *in  = src + c * esize;
        char       *out = (char *)comps[c] + row0 * esize;
        for (r = 0; r < rows; r++) {
            memcpy(out, in, esize);
            in  += ncomp * esize;
            out += esize;
        }
    }
}

/*
 * Shared validation for the two symmetric entry points. `comps` is inspected
 * as a const array either way; the callers cast their own qualification.
 */
static h5c_status_t check_interleaved(h5c_file_t *file, const char *path,
                                      const void *const *comps, size_t ncomp,
                                      size_t n, h5c_type_t type,
                                      size_t *elemsize, size_t dims[2])
{
    h5c_status_t st;
    size_t       i;

    if (ncomp == 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "ncomp is 0 for '%s'",
                         (path != NULL) ? path : "(null)");
    }
    dims[0] = n;
    dims[1] = ncomp;
    if ((st = h5c__check_common(file, path, 2, dims)) != H5C_OK) {
        return st;
    }
    if (comps == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "component pointer array is NULL for '%s'", path);
    }
    for (i = 0; i < ncomp; i++) {
        if (comps[i] == NULL) {
            return h5c__fail(H5C_ERR_INVALID_ARG,
                             "component %lu is NULL for '%s'",
                             (unsigned long)i, path);
        }
    }
    *elemsize = h5c_type_size(type);
    if (*elemsize == 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no interleaved mapping for '%s'",
                         (int)type, path);
    }
    return h5c__ensure_init();
}

/*
 * Selects rows [row0, row0 + rows) of the whole file space, all `ncomp`
 * columns. Returns the file dataspace (caller closes) or H5I_INVALID_HID.
 */
static hid_t select_rows(hid_t did, size_t row0, size_t rows, size_t ncomp)
{
    hsize_t start[2], count[2];
    hid_t   fsid;

    fsid = H5Dget_space(did);
    if (fsid < 0) {
        return H5I_INVALID_HID;
    }
    start[0] = (hsize_t)row0;
    start[1] = 0;
    count[0] = (hsize_t)rows;
    count[1] = (hsize_t)ncomp;
    if (H5Sselect_hyperslab(fsid, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        H5Sclose(fsid);
        return H5I_INVALID_HID;
    }
    return fsid;
}

/* Contiguous memory dataspace of `rows` x `ncomp`. */
static hid_t tile_memspace(size_t rows, size_t ncomp)
{
    hsize_t mdims[2];

    mdims[0] = (hsize_t)rows;
    mdims[1] = (hsize_t)ncomp;
    return H5Screate_simple(2, mdims, NULL);
}

/* ------------------------------------------------------------------ */
/* write                                                               */
/* ------------------------------------------------------------------ */

/* Opens or creates the [n, ncomp] dataset, honouring H5C_WRITE_REPLACE. */
static h5c_status_t open_for_write(h5c_file_t *file, const char *path,
                                   h5c_type_t type, const size_t dims[2],
                                   unsigned flags, hid_t *out_did)
{
    h5c_status_t st;
    hid_t        ftype, sid, did;
    int          existed;

    *out_did = H5I_INVALID_HID;

    ftype = h5c__file_type(type);
    if (ftype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping for '%s'",
                         (int)type, path);
    }

    existed = h5c_exists(file, path);
    if (existed && (flags & H5C_WRITE_REPLACE)) {
        if (H5Ldelete(file->fid, path, H5P_DEFAULT) < 0) {
            return h5c__fail_hdf5(-1, "cannot replace existing '%s'", path);
        }
        existed = 0;
    }

    if (existed) {
        h5c_dataset_info_t info;

        did = H5Dopen2(file->fid, path, H5P_DEFAULT);
        if (did < 0) {
            return h5c__fail_hdf5((long)did,
                                  "'%s' exists but is not a dataset", path);
        }
        if ((st = h5c__info_from_dset(did, &info)) != H5C_OK) {
            H5Dclose(did);
            return st;
        }
        if (h5c__shape_equals(&info, path, 2, dims) != H5C_OK) {
            H5Dclose(did);
            return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                             "'%s' already exists with a different shape; "
                             "pass H5C_WRITE_REPLACE to overwrite it", path);
        }
    } else {
        sid = h5c__make_space(2, dims);
        if (sid < 0) {
            return h5c__fail_hdf5((long)sid,
                                  "cannot build dataspace for '%s'", path);
        }
        did = H5Dcreate2(file->fid, path, ftype, sid,
                         h5c__lcpl(), H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(sid);
        if (did < 0) {
            return h5c__fail_hdf5((long)did, "cannot create '%s'", path);
        }
    }

    *out_did = did;
    return H5C_OK;
}

static h5c_status_t write_interleaved_impl(h5c_file_t *file, const char *path,
                                           const void *const *comps,
                                           size_t ncomp, size_t n,
                                           h5c_type_t type, unsigned flags)
{
    h5c_status_t st;
    size_t       dims[2], esize, row_bytes, rows, row0;
    hid_t        mtype, did = H5I_INVALID_HID;
    char        *stage;

    if ((st = check_interleaved(file, path, comps, ncomp, n, type,
                                &esize, dims)) != H5C_OK) {
        return st;
    }
    if (file->readonly) {
        return h5c__fail(H5C_ERR_STATE,
                         "file is open read-only, cannot write '%s'", path);
    }

    mtype = h5c__mem_type(type);
    if (mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping for '%s'",
                         (int)type, path);
    }
    if ((st = open_for_write(file, path, type, dims, flags, &did)) != H5C_OK) {
        return st;
    }

    row_bytes = ncomp * esize;
    rows      = h5c__tile_rows(n, row_bytes);

    stage = (char *)malloc(rows * row_bytes);
    if (stage == NULL) {
        H5Dclose(did);
        return h5c__fail(H5C_ERR_NOMEM,
                         "cannot allocate %lu staging bytes for '%s'",
                         (unsigned long)(rows * row_bytes), path);
    }

    st = H5C_OK;
    for (row0 = 0; row0 < n; row0 += rows) {
        size_t this_rows = (n - row0 < rows) ? (n - row0) : rows;
        hid_t  fsid, msid;

        h5c__pack_tile(stage, comps, ncomp, row0, this_rows, esize);

        fsid = select_rows(did, row0, this_rows, ncomp);
        if (fsid < 0) {
            st = h5c__fail_hdf5(-1, "cannot select rows of '%s'", path);
            break;
        }
        msid = tile_memspace(this_rows, ncomp);
        if (msid < 0) {
            H5Sclose(fsid);
            st = h5c__fail_hdf5(-1, "cannot build memory space for '%s'", path);
            break;
        }
        if (H5Dwrite(did, mtype, msid, fsid, H5P_DEFAULT, stage) < 0) {
            st = h5c__fail_hdf5(-1, "H5Dwrite failed for '%s'", path);
        }
        H5Sclose(msid);
        H5Sclose(fsid);
        if (st != H5C_OK) {
            break;
        }
    }

    free(stage);
    H5Dclose(did);
    return st;
}

h5c_status_t h5c_write_interleaved(h5c_file_t *file, const char *path,
                                   const void *const *comps, size_t ncomp,
                                   size_t n, h5c_type_t type, unsigned flags)
{
    return h5c__record(file, write_interleaved_impl(file, path, comps, ncomp,
                                                    n, type, flags));
}

/* ------------------------------------------------------------------ */
/* read                                                                */
/* ------------------------------------------------------------------ */

static h5c_status_t read_interleaved_impl(h5c_file_t *file, const char *path,
                                          void *const *comps, size_t ncomp,
                                          size_t n, h5c_type_t type)
{
    h5c_status_t       st;
    h5c_dataset_info_t info;
    size_t             dims[2], esize, row_bytes, rows, row0;
    hid_t              mtype, did;
    char              *stage;

    if ((st = check_interleaved(file, path, (const void *const *)comps,
                                ncomp, n, type, &esize, dims)) != H5C_OK) {
        return st;
    }

    mtype = h5c__mem_type(type);
    if (mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping for '%s'",
                         (int)type, path);
    }

    did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND, "no dataset at '%s'", path);
    }
    if ((st = h5c__info_from_dset(did, &info)) != H5C_OK) {
        H5Dclose(did);
        return st;
    }
    if ((st = h5c__shape_equals(&info, path, 2, dims)) != H5C_OK) {
        H5Dclose(did);
        return st;
    }

    row_bytes = ncomp * esize;
    rows      = h5c__tile_rows(n, row_bytes);

    stage = (char *)malloc(rows * row_bytes);
    if (stage == NULL) {
        H5Dclose(did);
        return h5c__fail(H5C_ERR_NOMEM,
                         "cannot allocate %lu staging bytes for '%s'",
                         (unsigned long)(rows * row_bytes), path);
    }

    st = H5C_OK;
    for (row0 = 0; row0 < n; row0 += rows) {
        size_t this_rows = (n - row0 < rows) ? (n - row0) : rows;
        hid_t  fsid, msid;

        fsid = select_rows(did, row0, this_rows, ncomp);
        if (fsid < 0) {
            st = h5c__fail_hdf5(-1, "cannot select rows of '%s'", path);
            break;
        }
        msid = tile_memspace(this_rows, ncomp);
        if (msid < 0) {
            H5Sclose(fsid);
            st = h5c__fail_hdf5(-1, "cannot build memory space for '%s'", path);
            break;
        }
        if (H5Dread(did, mtype, msid, fsid, H5P_DEFAULT, stage) < 0) {
            st = h5c__fail_hdf5(-1, "H5Dread failed for '%s'", path);
        }
        H5Sclose(msid);
        H5Sclose(fsid);
        if (st != H5C_OK) {
            break;
        }
        h5c__unpack_tile(stage, comps, ncomp, row0, this_rows, esize);
    }

    free(stage);
    H5Dclose(did);
    return st;
}

h5c_status_t h5c_read_interleaved(h5c_file_t *file, const char *path,
                                  void *const *comps, size_t ncomp,
                                  size_t n, h5c_type_t type)
{
    return h5c__record(file, read_interleaved_impl(file, path, comps, ncomp,
                                                   n, type));
}

/* ------------------------------------------------------------------ */
/* single component                                                    */
/* ------------------------------------------------------------------ */

static h5c_status_t read_component_impl(h5c_file_t *file, const char *path,
                                        void *buf, size_t comp, size_t n,
                                        h5c_type_t type)
{
    h5c_status_t       st;
    h5c_dataset_info_t info;
    hsize_t            start[2], count[2], mdim;
    size_t             dims1;
    hid_t              mtype, did, fsid, msid;

    dims1 = n;
    if ((st = h5c__check_common(file, path, 1, &dims1)) != H5C_OK) {
        return st;
    }
    if (buf == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "buffer is NULL for '%s'", path);
    }
    if (h5c_type_size(type) == 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no interleaved mapping for '%s'",
                         (int)type, path);
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    mtype = h5c__mem_type(type);
    if (mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping for '%s'",
                         (int)type, path);
    }

    did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND, "no dataset at '%s'", path);
    }
    if ((st = h5c__info_from_dset(did, &info)) != H5C_OK) {
        H5Dclose(did);
        return st;
    }
    /* ncomp comes from the stored dataset, never from the caller. */
    if (info.rank != 2) {
        H5Dclose(did);
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         "'%s' has rank %d, expected an [n, ncomp] dataset",
                         path, info.rank);
    }
    if (info.dims[0] != n) {
        H5Dclose(did);
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         "'%s' has %lu rows, expected %lu",
                         path, (unsigned long)info.dims[0], (unsigned long)n);
    }
    if (comp >= info.dims[1]) {
        H5Dclose(did);
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "component %lu is out of range for '%s' (ncomp %lu)",
                         (unsigned long)comp, path,
                         (unsigned long)info.dims[1]);
    }

    /*
     * One column of the file: start = {0, comp}, count = {n, 1}. Strided in
     * the file, contiguous in memory, and only 1/ncomp of the bytes move.
     */
    fsid = H5Dget_space(did);
    if (fsid < 0) {
        H5Dclose(did);
        return h5c__fail_hdf5((long)fsid, "H5Dget_space failed for '%s'", path);
    }
    start[0] = 0;
    start[1] = (hsize_t)comp;
    count[0] = (hsize_t)n;
    count[1] = 1;
    if (H5Sselect_hyperslab(fsid, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        H5Sclose(fsid);
        H5Dclose(did);
        return h5c__fail_hdf5(-1, "cannot select component %lu of '%s'",
                              (unsigned long)comp, path);
    }

    mdim = (hsize_t)n;
    msid = H5Screate_simple(1, &mdim, NULL);
    if (msid < 0) {
        H5Sclose(fsid);
        H5Dclose(did);
        return h5c__fail_hdf5((long)msid,
                              "cannot build memory space for '%s'", path);
    }

    st = H5C_OK;
    if (H5Dread(did, mtype, msid, fsid, H5P_DEFAULT, buf) < 0) {
        st = h5c__fail_hdf5(-1, "H5Dread failed for component %lu of '%s'",
                            (unsigned long)comp, path);
    }
    H5Sclose(msid);
    H5Sclose(fsid);
    H5Dclose(did);
    return st;
}

h5c_status_t h5c_read_component(h5c_file_t *file, const char *path,
                                void *buf, size_t comp, size_t n,
                                h5c_type_t type)
{
    return h5c__record(file, read_component_impl(file, path, buf, comp, n, type));
}
