#include "h5c_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

h5c_status_t h5c__check_common(h5c_file_t *file, const char *path,
                                 int rank, const size_t *dims)
{
    int i;

    if (file == NULL || file->fid < 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "file handle is NULL or closed");
    }
    if (path == NULL || path[0] == '\0') {
        return h5c__fail(H5C_ERR_INVALID_ARG, "empty dataset path");
    }
    if (rank < 0 || rank > H5C_MAX_RANK) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "rank %d outside [0, %d] for '%s'",
                         rank, H5C_MAX_RANK, path);
    }
    if (rank > 0 && dims == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "dims is NULL for '%s'", path);
    }
    /*
     * A zero extent is deliberately allowed. Parallel I/O needs it (a rank may
     * own no rows) and a zero-extent dataset is a valid HDF5 object, so the
     * serial and parallel paths share one policy instead of two. Callers get
     * no "uninitialised n" guard here; that trade is intentional.
     */
    (void)i;
    return H5C_OK;
}

size_t h5c__count(int rank, const size_t *dims)
{
    size_t n = 1;
    int    i;

    for (i = 0; i < rank; i++) {
        n *= dims[i];
    }
    return n;
}

/* Creates a dataspace for `rank`/`dims`. rank 0 gives a scalar space. */
hid_t h5c__make_space(int rank, const size_t *dims)
{
    hsize_t hdims[H5C_MAX_RANK];
    int     i;

    if (rank == 0) {
        return H5Screate(H5S_SCALAR);
    }
    for (i = 0; i < rank; i++) {
        hdims[i] = (hsize_t)dims[i];
    }
    return H5Screate_simple(rank, hdims, NULL);
}

/* Fills `out` from an already-open dataset. */
h5c_status_t h5c__info_from_dset(hid_t did, h5c_dataset_info_t *out)
{
    hsize_t hdims[H5C_MAX_RANK];
    hid_t   sid, tid;
    int     ndims, i;

    memset(out, 0, sizeof *out);
    out->count = 1;

    sid = H5Dget_space(did);
    if (sid < 0) {
        return h5c__fail_hdf5((long)sid, "H5Dget_space failed");
    }
    ndims = H5Sget_simple_extent_ndims(sid);
    if (ndims < 0) {
        H5Sclose(sid);
        return h5c__fail_hdf5((long)ndims, "H5Sget_simple_extent_ndims failed");
    }
    if (ndims > H5C_MAX_RANK) {
        H5Sclose(sid);
        return h5c__fail(H5C_ERR_UNSUPPORTED,
                         "dataset rank %d exceeds H5C_MAX_RANK", ndims);
    }
    if (ndims > 0 && H5Sget_simple_extent_dims(sid, hdims, NULL) < 0) {
        H5Sclose(sid);
        return h5c__fail_hdf5(-1, "H5Sget_simple_extent_dims failed");
    }
    H5Sclose(sid);

    out->rank = ndims;
    for (i = 0; i < ndims; i++) {
        out->dims[i] = (size_t)hdims[i];
        out->count  *= (size_t)hdims[i];
    }

    tid = H5Dget_type(did);
    if (tid < 0) {
        return h5c__fail_hdf5((long)tid, "H5Dget_type failed");
    }
    out->type = h5c__type_from_hid(tid);
    H5Tclose(tid);

    return H5C_OK;
}

/* Compares a stored shape against the caller's expectation. */
h5c_status_t h5c__shape_equals(const h5c_dataset_info_t *info,
                                 const char *path, int rank, const size_t *dims)
{
    int i;

    if (info->rank != rank) {
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         "'%s' has rank %d, expected %d",
                         path, info->rank, rank);
    }
    for (i = 0; i < rank; i++) {
        if (info->dims[i] != dims[i]) {
            return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                             "'%s' dims[%d] is %lu, expected %lu",
                             path, i,
                             (unsigned long)info->dims[i],
                             (unsigned long)dims[i]);
        }
    }
    return H5C_OK;
}

int h5c_exists(h5c_file_t *file, const char *path)
{
    htri_t r;

    if (file == NULL || file->fid < 0 || path == NULL || path[0] == '\0') {
        return 0;
    }
    r = H5Lexists(file->fid, path, H5P_DEFAULT);
    return (r > 0) ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* write                                                               */
/* ------------------------------------------------------------------ */

static h5c_status_t write_impl(h5c_file_t *file, const char *path,
                               const void *buf, h5c_type_t type,
                               int rank, const size_t *dims, unsigned flags)
{
    h5c_status_t st;
    hid_t        ftype, mtype, sid = H5I_INVALID_HID, did = H5I_INVALID_HID;
    int          existed;

    if ((st = h5c__check_common(file, path, rank, dims)) != H5C_OK) {
        return st;
    }
    if (buf == NULL && h5c__count(rank, dims) > 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "buffer is NULL for '%s'", path);
    }
    if (file->readonly) {
        return h5c__fail(H5C_ERR_STATE,
                         "file is open read-only, cannot write '%s'", path);
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    ftype = h5c__file_type(type);
    mtype = h5c__mem_type(type);
    if (ftype == H5I_INVALID_HID || mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping (use the string API)",
                         (int)type);
    }

    existed = h5c_exists(file, path);

    if (existed && (flags & H5C_WRITE_REPLACE)) {
        if (H5Ldelete(file->fid, path, H5P_DEFAULT) < 0) {
            return h5c__fail_hdf5(-1, "cannot replace existing '%s'", path);
        }
        existed = 0;
    }

    if (existed) {
        /* Write in place; the stored shape must match exactly. */
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
        if ((st = h5c__shape_equals(&info, path, rank, dims)) != H5C_OK) {
            H5Dclose(did);
            return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                             "'%s' already exists with a different shape; "
                             "pass H5C_WRITE_REPLACE to overwrite it", path);
        }
    } else {
        sid = h5c__make_space(rank, dims);
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

    /* An empty dataset exists but carries no elements: nothing to transfer. */
    if (h5c__count(rank, dims) > 0 &&
        H5Dwrite(did, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) {
        H5Dclose(did);
        return h5c__fail_hdf5(-1, "H5Dwrite failed for '%s'", path);
    }
    H5Dclose(did);
    return H5C_OK;
}

h5c_status_t h5c_write(h5c_file_t *file, const char *path, const void *buf,
                       h5c_type_t type, int rank, const size_t *dims,
                       unsigned flags)
{
    return h5c__record(file, write_impl(file, path, buf, type, rank, dims, flags));
}

/* ------------------------------------------------------------------ */
/* read                                                                */
/* ------------------------------------------------------------------ */

static h5c_status_t read_impl(h5c_file_t *file, const char *path, void *buf,
                              h5c_type_t type, int rank, const size_t *dims)
{
    h5c_status_t       st;
    h5c_dataset_info_t info;
    hid_t              mtype, did;

    if ((st = h5c__check_common(file, path, rank, dims)) != H5C_OK) {
        return st;
    }
    if (buf == NULL && h5c__count(rank, dims) > 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "buffer is NULL for '%s'", path);
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    mtype = h5c__mem_type(type);
    if (mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping (use the string API)",
                         (int)type);
    }

    did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND, "no dataset at '%s'", path);
    }
    if ((st = h5c__info_from_dset(did, &info)) != H5C_OK) {
        H5Dclose(did);
        return st;
    }
    if ((st = h5c__shape_equals(&info, path, rank, dims)) != H5C_OK) {
        H5Dclose(did);
        return st;
    }
    /* HDF5 converts between the stored type and `mtype` where it can. */
    if (info.count > 0 &&
        H5Dread(did, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) {
        H5Dclose(did);
        return h5c__fail_hdf5(-1, "H5Dread failed for '%s'", path);
    }
    H5Dclose(did);
    return H5C_OK;
}

h5c_status_t h5c_read(h5c_file_t *file, const char *path, void *buf,
                      h5c_type_t type, int rank, const size_t *dims)
{
    return h5c__record(file, read_impl(file, path, buf, type, rank, dims));
}

/* ------------------------------------------------------------------ */
/* metadata and allocating read                                        */
/* ------------------------------------------------------------------ */

static h5c_status_t info_impl(h5c_file_t *file, const char *path,
                              h5c_dataset_info_t *out)
{
    h5c_status_t st;
    hid_t        did;

    if (file == NULL || file->fid < 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "file handle is NULL or closed");
    }
    if (path == NULL || path[0] == '\0' || out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_dataset_info: bad argument");
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND, "no dataset at '%s'", path);
    }
    st = h5c__info_from_dset(did, out);
    H5Dclose(did);
    return st;
}

h5c_status_t h5c_dataset_info(h5c_file_t *file, const char *path,
                              h5c_dataset_info_t *out)
{
    return h5c__record(file, info_impl(file, path, out));
}

static h5c_status_t read_alloc_impl(h5c_file_t *file, const char *path,
                                    h5c_type_t type, void **buf,
                                    h5c_dataset_info_t *info)
{
    h5c_status_t       st;
    h5c_dataset_info_t local;
    size_t             esize;
    void              *mem;

    if (buf == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_read_alloc: buf is NULL");
    }
    *buf = NULL;

    esize = h5c_type_size(type);
    if (esize == 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "h5c_read_alloc does not support type %d", (int)type);
    }
    if ((st = info_impl(file, path, &local)) != H5C_OK) {
        return st;
    }

    /* Allocate at least one byte so *buf == NULL always means failure. */
    mem = malloc(local.count > 0 ? local.count * esize : 1);
    if (mem == NULL) {
        return h5c__fail(H5C_ERR_NOMEM,
                         "cannot allocate %lu bytes for '%s'",
                         (unsigned long)(local.count * esize), path);
    }

    st = read_impl(file, path, mem, type, local.rank,
                   (local.rank > 0) ? local.dims : NULL);
    if (st != H5C_OK) {
        free(mem);
        return st;
    }

    *buf = mem;
    if (info != NULL) {
        *info = local;
    }
    return H5C_OK;
}

h5c_status_t h5c_read_alloc(h5c_file_t *file, const char *path,
                            h5c_type_t type, void **buf,
                            h5c_dataset_info_t *info)
{
    return h5c__record(file, read_alloc_impl(file, path, type, buf, info));
}
