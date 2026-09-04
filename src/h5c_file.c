#include "h5c_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Whether `path` can be opened for reading at all. Deliberately plain stdio
 * rather than access(): the library is built as strict C11, where POSIX
 * declarations are not visible without a feature macro, and this only needs
 * to separate "absent" from "present but unusable".
 */
static int path_exists(const char *path)
{
    FILE *fp = fopen(path, "rb");

    if (fp == NULL) {
        return 0;
    }
    fclose(fp);
    return 1;
}

h5c_status_t h5c_open(const char *path, h5c_mode_t mode, h5c_file_t **out)
{
    h5c_status_t st;
    h5c_file_t  *file;
    hid_t        fid;

    if (out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_open: out is NULL");
    }
    *out = NULL;
    if (path == NULL || path[0] == '\0') {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_open: empty path");
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    switch (mode) {
    case H5C_TRUNCATE:
        fid = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        break;
    case H5C_READ:
        fid = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
        break;
    case H5C_READWRITE:
        fid = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
        break;
    default:
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_open: bad mode %d", (int)mode);
    }

    if (fid < 0) {
        /*
         * Distinguish "not there" from "there but unusable". h5c.h documents
         * H5C_ERR_NOT_FOUND for a missing file, and a caller deciding whether
         * to create a fresh file needs that apart from a corrupt or
         * permission-denied one, which stays H5C_ERR_HDF5.
         */
        if (mode != H5C_TRUNCATE && !path_exists(path)) {
            return h5c__fail(H5C_ERR_NOT_FOUND, "no such file: '%s'", path);
        }
        return h5c__fail_hdf5((long)fid, "cannot open '%s' (mode %d)",
                              path, (int)mode);
    }

    file = (h5c_file_t *)calloc(1, sizeof *file);
    if (file == NULL) {
        H5Fclose(fid);
        return h5c__fail(H5C_ERR_NOMEM, "h5c_open: allocation failed");
    }
    file->fid      = fid;
    file->sticky   = H5C_OK;
    file->borrowed = 0;
    file->readonly = (mode == H5C_READ);

    *out = file;
    return H5C_OK;
}

h5c_status_t h5c_file_from_hid(hid_t fid, h5c_file_t **out)
{
    h5c_file_t *file;

    if (out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_file_from_hid: out is NULL");
    }
    *out = NULL;
    if (fid < 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_file_from_hid: invalid id");
    }

    file = (h5c_file_t *)calloc(1, sizeof *file);
    if (file == NULL) {
        return h5c__fail(H5C_ERR_NOMEM, "h5c_file_from_hid: allocation failed");
    }
    file->fid      = fid;
    file->sticky   = H5C_OK;
    file->borrowed = 1;
    file->readonly = 0;

    *out = file;
    return H5C_OK;
}

h5c_status_t h5c_close(h5c_file_t *file)
{
    h5c_status_t st = H5C_OK;

    if (file == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_close: file is NULL");
    }

    /*
     * Reports ONLY whether closing succeeded. It deliberately does not report
     * the sticky status: a caller must be able to tell "the file may not have
     * been flushed" from "some earlier, already-handled call failed". Query
     * h5c_file_status() before closing if you use the sticky pattern.
     */
    if (!file->borrowed && file->fid >= 0) {
        if (H5Fclose(file->fid) < 0) {
            st = h5c__fail_hdf5(-1, "H5Fclose failed");
        }
    }
    free(file);
    return st;
}

h5c_status_t h5c_file_status(const h5c_file_t *file)
{
    return (file == NULL) ? H5C_ERR_INVALID_ARG : file->sticky;
}

void h5c_file_clear_status(h5c_file_t *file)
{
    if (file != NULL) {
        file->sticky = H5C_OK;
    }
}

int h5c_is_parallel(const h5c_file_t *file)
{
    /*
     * Deliberately lives here rather than in h5c_parallel.c: h5c.h must stay
     * usable without <mpi.h>, and a serial build has no parallel state at all.
     */
#ifdef H5C_HAVE_PARALLEL
    return (file != NULL && file->parallel) ? 1 : 0;
#else
    (void)file;
    return 0;
#endif
}

hid_t h5c_file_hid(const h5c_file_t *file)
{
    return (file == NULL) ? H5I_INVALID_HID : file->fid;
}

void h5c_free(void *ptr)
{
    free(ptr);
}
