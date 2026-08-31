/*
 * Scalar string datasets.
 *
 * The stored representation deliberately matches what h5fortran produces:
 * a SCALAR dataset of fixed-length H5T_C_S1 with STRSIZE == strlen(value)
 * (at least 1), STRPAD H5T_STR_SPACEPAD and CSET ASCII. h5fortran writes
 * `trim(str)` into a type sized by len_trim(), so the padding character on
 * the wire is a space, never a NUL.
 *
 * Reading accepts both fixed-length and variable-length strings and always
 * hands back a NUL-terminated buffer. Trailing spaces and NULs of a
 * fixed-length value are stripped, which is what SPACEPAD means.
 */
#include "h5c_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* shared helpers (also used by h5c_attribute.c)                       */
/* ------------------------------------------------------------------ */

/*
 * Declared here rather than in h5c_internal.h because that header is shared
 * with other work in progress. Keep these in sync with h5c_attribute.c.
 */
hid_t        h5c__str_fixed_type(size_t len);
char        *h5c__str_pad_buffer(const char *value, size_t *len_out);
h5c_status_t h5c__str_read_id(hid_t id, int is_attr, const char *what,
                              char **out);

/* Fixed-length H5T_C_S1 with SPACEPAD. The caller closes the returned id. */
hid_t h5c__str_fixed_type(size_t len)
{
    hid_t tid;

    if (len == 0) {
        len = 1; /* HDF5 rejects a zero-sized string type */
    }
    tid = H5Tcopy(H5T_C_S1);
    if (tid < 0) {
        h5c__fail_hdf5((long)tid, "H5Tcopy(H5T_C_S1) failed");
        return H5I_INVALID_HID;
    }
    if (H5Tset_size(tid, len) < 0 ||
        H5Tset_strpad(tid, H5T_STR_SPACEPAD) < 0 ||
        H5Tset_cset(tid, H5T_CSET_ASCII) < 0) {
        H5Tclose(tid);
        h5c__fail_hdf5(-1, "cannot configure a fixed-length string type");
        return H5I_INVALID_HID;
    }
    return tid;
}

/*
 * Staging buffer of exactly *len_out bytes, space padded and NOT
 * NUL-terminated, as the fixed-length type expects. Released with free().
 */
char *h5c__str_pad_buffer(const char *value, size_t *len_out)
{
    size_t n, len;
    char  *buf;

    len = strlen(value);
    n   = (len > 0) ? len : 1; /* the empty string becomes a single space */

    buf = (char *)malloc(n);
    if (buf == NULL) {
        return NULL;
    }
    memset(buf, ' ', n);
    memcpy(buf, value, len);

    *len_out = n;
    return buf;
}

/* Trims trailing spaces and NULs in place. */
static void trim_right(char *s, size_t n)
{
    while (n > 0 && (s[n - 1] == ' ' || s[n - 1] == '\0')) {
        n--;
    }
    s[n] = '\0';
}

/*
 * Reads the scalar string held by an attribute (`is_attr`) or a dataset id.
 * `*out` is a malloc'd NUL-terminated buffer on success and NULL otherwise.
 */
h5c_status_t h5c__str_read_id(hid_t id, int is_attr, const char *what,
                              char **out)
{
    h5c_status_t st    = H5C_OK;
    hid_t        ftype = H5I_INVALID_HID;
    hid_t        mtype = H5I_INVALID_HID;
    hid_t        sid   = H5I_INVALID_HID;
    char        *buf   = NULL;
    size_t       n;

    *out  = NULL;
    ftype = is_attr ? H5Aget_type(id) : H5Dget_type(id);
    if (ftype < 0) {
        return h5c__fail_hdf5((long)ftype, "cannot query the type of '%s'",
                              what);
    }
    if (H5Tget_class(ftype) != H5T_STRING) {
        H5Tclose(ftype);
        return h5c__fail(H5C_ERR_TYPE_MISMATCH, "'%s' is not a string", what);
    }

    sid = is_attr ? H5Aget_space(id) : H5Dget_space(id);
    if (sid < 0) {
        H5Tclose(ftype);
        return h5c__fail_hdf5((long)sid, "cannot query the space of '%s'",
                              what);
    }
    if (H5Sget_simple_extent_npoints(sid) != 1) {
        H5Sclose(sid);
        H5Tclose(ftype);
        return h5c__fail(H5C_ERR_UNSUPPORTED,
                         "'%s' holds more than one string", what);
    }

    if (H5Tis_variable_str(ftype) > 0) {
        char *raw = NULL;

        mtype = H5Tcopy(H5T_C_S1);
        if (mtype < 0 || H5Tset_size(mtype, H5T_VARIABLE) < 0 ||
            H5Tset_cset(mtype, H5Tget_cset(ftype)) < 0) {
            st = h5c__fail_hdf5(-1,
                                "cannot build a variable-length string type");
            goto done;
        }
        if ((is_attr ? H5Aread(id, mtype, &raw)
                     : H5Dread(id, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               &raw)) < 0) {
            st = h5c__fail_hdf5(-1, "cannot read '%s'", what);
            goto done;
        }

        n   = (raw != NULL) ? strlen(raw) : 0;
        buf = (char *)malloc(n + 1);
        if (buf != NULL) {
            memcpy(buf, (raw != NULL) ? raw : "", n);
            buf[n] = '\0';
        }
        H5Treclaim(mtype, sid, H5P_DEFAULT, &raw);
        if (buf == NULL) {
            st = h5c__fail(H5C_ERR_NOMEM, "cannot allocate %lu bytes for '%s'",
                           (unsigned long)(n + 1), what);
            goto done;
        }
        /* A variable-length value is exact: there is no padding to strip. */
    } else {
        n = H5Tget_size(ftype);
        if (n == 0) {
            st = h5c__fail_hdf5(-1, "H5Tget_size failed for '%s'", what);
            goto done;
        }
        /* Read through a copy of the stored type: no conversion, no surprises. */
        mtype = H5Tcopy(ftype);
        if (mtype < 0) {
            st = h5c__fail_hdf5((long)mtype, "H5Tcopy failed for '%s'", what);
            goto done;
        }
        buf = (char *)malloc(n + 1);
        if (buf == NULL) {
            st = h5c__fail(H5C_ERR_NOMEM, "cannot allocate %lu bytes for '%s'",
                           (unsigned long)(n + 1), what);
            goto done;
        }
        if ((is_attr ? H5Aread(id, mtype, buf)
                     : H5Dread(id, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               buf)) < 0) {
            free(buf);
            buf = NULL;
            st  = h5c__fail_hdf5(-1, "cannot read '%s'", what);
            goto done;
        }
        buf[n] = '\0';
        trim_right(buf, n);
    }

    *out = buf;
    buf  = NULL;

done:
    if (buf != NULL) {
        free(buf);
    }
    if (mtype >= 0) {
        H5Tclose(mtype);
    }
    if (sid >= 0) {
        H5Sclose(sid);
    }
    if (ftype >= 0) {
        H5Tclose(ftype);
    }
    return st;
}

/* ------------------------------------------------------------------ */
/* write                                                               */
/* ------------------------------------------------------------------ */

static h5c_status_t check_write(h5c_file_t *file, const char *path,
                                const char *value)
{
    h5c_status_t st;

    if ((st = h5c__check_common(file, path, 0, NULL)) != H5C_OK) {
        return st;
    }
    if (value == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "value is NULL for '%s'", path);
    }
    if (file->readonly) {
        return h5c__fail(H5C_ERR_STATE,
                         "file is open read-only, cannot write '%s'", path);
    }
    return h5c__ensure_init();
}

/*
 * Removes an existing link when replacing. Otherwise reports whether the
 * stored string can still be written in place, which requires exactly the
 * same length: a fixed-length datatype cannot grow.
 */
static h5c_status_t prepare_target(h5c_file_t *file, const char *path,
                                   size_t len, unsigned flags, int *existed)
{
    hid_t  did, tid;
    size_t stored;

    *existed = h5c_exists(file, path);
    if (!*existed) {
        return H5C_OK;
    }
    if (flags & H5C_WRITE_REPLACE) {
        if (H5Ldelete(file->fid, path, H5P_DEFAULT) < 0) {
            return h5c__fail_hdf5(-1, "cannot replace existing '%s'", path);
        }
        *existed = 0;
        return H5C_OK;
    }

    did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail_hdf5((long)did, "'%s' exists but is not a dataset",
                              path);
    }
    tid = H5Dget_type(did);
    if (tid < 0) {
        H5Dclose(did);
        return h5c__fail_hdf5((long)tid, "H5Dget_type failed for '%s'", path);
    }
    stored = (H5Tget_class(tid) == H5T_STRING && H5Tis_variable_str(tid) <= 0)
                 ? H5Tget_size(tid)
                 : 0;
    H5Tclose(tid);
    H5Dclose(did);

    if (stored != len) {
        return h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                         "'%s' already stores a string of a different length; "
                         "pass H5C_WRITE_REPLACE to overwrite it", path);
    }
    return H5C_OK;
}

static h5c_status_t write_fixed_impl(h5c_file_t *file, const char *path,
                                     const char *value, unsigned flags)
{
    h5c_status_t st;
    hid_t        tid = H5I_INVALID_HID, sid = H5I_INVALID_HID;
    hid_t        did = H5I_INVALID_HID;
    char        *buf = NULL;
    size_t       len = 0;
    int          existed = 0;

    if ((st = check_write(file, path, value)) != H5C_OK) {
        return st;
    }

    buf = h5c__str_pad_buffer(value, &len);
    if (buf == NULL) {
        return h5c__fail(H5C_ERR_NOMEM, "cannot stage the value for '%s'",
                         path);
    }
    if ((st = prepare_target(file, path, len, flags, &existed)) != H5C_OK) {
        free(buf);
        return st;
    }

    tid = h5c__str_fixed_type(len);
    if (tid == H5I_INVALID_HID) {
        free(buf);
        return H5C_ERR_HDF5;
    }

    if (existed) {
        did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    } else {
        sid = H5Screate(H5S_SCALAR);
        if (sid < 0) {
            free(buf);
            H5Tclose(tid);
            return h5c__fail_hdf5((long)sid,
                                  "cannot create a scalar space for '%s'",
                                  path);
        }
        did = H5Dcreate2(file->fid, path, tid, sid,
                         h5c__lcpl(), H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(sid);
    }
    if (did < 0) {
        free(buf);
        H5Tclose(tid);
        return h5c__fail_hdf5((long)did, "cannot create '%s'", path);
    }

    st = H5C_OK;
    if (H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) {
        st = h5c__fail_hdf5(-1, "H5Dwrite failed for '%s'", path);
    }

    H5Dclose(did);
    H5Tclose(tid);
    free(buf);
    return st;
}

static h5c_status_t write_vlen_impl(h5c_file_t *file, const char *path,
                                    const char *value, unsigned flags)
{
    h5c_status_t st;
    hid_t        tid = H5I_INVALID_HID, sid = H5I_INVALID_HID;
    hid_t        did = H5I_INVALID_HID;

    if ((st = check_write(file, path, value)) != H5C_OK) {
        return st;
    }
    if (h5c_exists(file, path)) {
        /* A variable-length value cannot be rewritten in place safely. */
        if (!(flags & H5C_WRITE_REPLACE)) {
            return h5c__fail(H5C_ERR_EXISTS,
                             "'%s' already exists; pass H5C_WRITE_REPLACE",
                             path);
        }
        if (H5Ldelete(file->fid, path, H5P_DEFAULT) < 0) {
            return h5c__fail_hdf5(-1, "cannot replace existing '%s'", path);
        }
    }

    tid = H5Tcopy(H5T_C_S1);
    if (tid < 0 || H5Tset_size(tid, H5T_VARIABLE) < 0 ||
        H5Tset_cset(tid, H5T_CSET_ASCII) < 0) {
        if (tid >= 0) {
            H5Tclose(tid);
        }
        return h5c__fail_hdf5(-1, "cannot build a variable-length string type");
    }

    sid = H5Screate(H5S_SCALAR);
    if (sid < 0) {
        H5Tclose(tid);
        return h5c__fail_hdf5((long)sid,
                              "cannot create a scalar space for '%s'", path);
    }
    did = H5Dcreate2(file->fid, path, tid, sid,
                     h5c__lcpl(), H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(sid);
    if (did < 0) {
        H5Tclose(tid);
        return h5c__fail_hdf5((long)did, "cannot create '%s'", path);
    }

    st = H5C_OK;
    if (H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0) {
        st = h5c__fail_hdf5(-1, "H5Dwrite failed for '%s'", path);
    }

    H5Dclose(did);
    H5Tclose(tid);
    return st;
}

h5c_status_t h5c_write_string(h5c_file_t *file, const char *path,
                              const char *value, unsigned flags)
{
    return h5c__record(file, write_fixed_impl(file, path, value, flags));
}

h5c_status_t h5c_write_string_vlen(h5c_file_t *file, const char *path,
                                   const char *value, unsigned flags)
{
    return h5c__record(file, write_vlen_impl(file, path, value, flags));
}

/* ------------------------------------------------------------------ */
/* read                                                                */
/* ------------------------------------------------------------------ */

static h5c_status_t read_string_impl(h5c_file_t *file, const char *path,
                                     char **out)
{
    h5c_status_t st;
    hid_t        did;

    if (out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_read_string: out is NULL");
    }
    *out = NULL;

    if ((st = h5c__check_common(file, path, 0, NULL)) != H5C_OK) {
        return st;
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    did = H5Dopen2(file->fid, path, H5P_DEFAULT);
    if (did < 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND, "no dataset at '%s'", path);
    }
    st = h5c__str_read_id(did, 0, path, out);
    H5Dclose(did);
    return st;
}

h5c_status_t h5c_read_string(h5c_file_t *file, const char *path, char **out)
{
    return h5c__record(file, read_string_impl(file, path, out));
}

void h5c_free_string(char *s)
{
    free(s);
}
