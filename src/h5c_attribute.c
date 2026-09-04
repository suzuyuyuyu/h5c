/*
 * Attributes on datasets, groups and the root group.
 *
 * `obj_path` is resolved with H5Oopen, so it may name any object, including
 * "/" for the root group. String attributes use the same fixed-length,
 * space-padded H5T_C_S1 representation as h5c_write_string(), which is what
 * h5fortran's h5fort_write_attribute writes.
 *
 * Writing an attribute that already exists replaces it (delete then create),
 * matching h5fortran. Array attributes are out of scope by design: large data
 * belongs in a dataset.
 */
#include "h5c_internal.h"

#include <stdlib.h>
#include <string.h>

/* Defined in h5c_string.c. Keep these declarations in sync with that file. */
hid_t        h5c__str_fixed_type(size_t len);
char        *h5c__str_pad_buffer(const char *value, size_t *len_out);
h5c_status_t h5c__str_read_id(hid_t id, int is_attr, const char *what,
                              char **out);

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static h5c_status_t check_args(h5c_file_t *file, const char *obj_path,
                               const char *name)
{
    if (file == NULL || file->fid < 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "file handle is NULL or closed");
    }
    if (obj_path == NULL || obj_path[0] == '\0') {
        return h5c__fail(H5C_ERR_INVALID_ARG, "empty object path");
    }
    if (name == NULL || name[0] == '\0') {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "empty attribute name for '%s'", obj_path);
    }
    return h5c__ensure_init();
}

/* Opens `obj_path`. The caller closes the id with H5Oclose. */
static h5c_status_t open_object(h5c_file_t *file, const char *obj_path,
                                hid_t *out)
{
    hid_t oid = H5Oopen(file->fid, obj_path, H5P_DEFAULT);

    if (oid < 0) {
        *out = H5I_INVALID_HID;
        return h5c__fail(H5C_ERR_NOT_FOUND, "no object at '%s'", obj_path);
    }
    *out = oid;
    return H5C_OK;
}

/* Deletes an existing attribute of the same name, as h5fortran does. */
static h5c_status_t drop_existing(hid_t oid, const char *obj_path,
                                  const char *name)
{
    htri_t present = H5Aexists(oid, name);

    if (present < 0) {
        return h5c__fail_hdf5((long)present,
                              "H5Aexists failed for '%s' on '%s'",
                              name, obj_path);
    }
    if (present > 0 && H5Adelete(oid, name) < 0) {
        return h5c__fail_hdf5(-1, "cannot replace attribute '%s' on '%s'",
                              name, obj_path);
    }
    return H5C_OK;
}

/* Opens an existing attribute for reading. The caller closes it. */
static h5c_status_t open_attr(hid_t oid, const char *obj_path,
                              const char *name, hid_t *out)
{
    htri_t present;
    hid_t  aid;

    *out    = H5I_INVALID_HID;
    present = H5Aexists(oid, name);
    if (present < 0) {
        return h5c__fail_hdf5((long)present,
                              "H5Aexists failed for '%s' on '%s'",
                              name, obj_path);
    }
    if (present == 0) {
        return h5c__fail(H5C_ERR_NOT_FOUND, "no attribute '%s' on '%s'",
                         name, obj_path);
    }
    aid = H5Aopen(oid, name, H5P_DEFAULT);
    if (aid < 0) {
        return h5c__fail_hdf5((long)aid, "cannot open attribute '%s' on '%s'",
                              name, obj_path);
    }
    *out = aid;
    return H5C_OK;
}

/* ------------------------------------------------------------------ */
/* string attributes                                                   */
/* ------------------------------------------------------------------ */

static h5c_status_t write_str_impl(h5c_file_t *file, const char *obj_path,
                                   const char *name, const char *value)
{
    h5c_status_t st;
    hid_t        oid = H5I_INVALID_HID, tid = H5I_INVALID_HID;
    hid_t        sid = H5I_INVALID_HID, aid = H5I_INVALID_HID;
    char        *buf = NULL;
    size_t       len = 0;

    if ((st = check_args(file, obj_path, name)) != H5C_OK) {
        return st;
    }
    if (value == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "value is NULL for attribute '%s'", name);
    }
    if (file->readonly) {
        return h5c__fail(H5C_ERR_STATE,
                         "file is open read-only, cannot write '%s'", name);
    }
    if ((st = open_object(file, obj_path, &oid)) != H5C_OK) {
        return st;
    }

    buf = h5c__str_pad_buffer(value, &len);
    if (buf == NULL) {
        H5Oclose(oid);
        return h5c__fail(H5C_ERR_NOMEM,
                         "cannot stage the value of attribute '%s'", name);
    }

    tid = h5c__str_fixed_type(len);
    if (tid == H5I_INVALID_HID) {
        free(buf);
        H5Oclose(oid);
        return H5C_ERR_HDF5;
    }
    sid = H5Screate(H5S_SCALAR);
    if (sid < 0) {
        st = h5c__fail_hdf5((long)sid,
                            "cannot create a scalar space for '%s'", name);
        goto done;
    }
    if ((st = drop_existing(oid, obj_path, name)) != H5C_OK) {
        goto done;
    }

    aid = H5Acreate2(oid, name, tid, sid, H5P_DEFAULT, H5P_DEFAULT);
    if (aid < 0) {
        st = h5c__fail_hdf5((long)aid, "cannot create attribute '%s' on '%s'",
                            name, obj_path);
        goto done;
    }
    if (H5Awrite(aid, tid, buf) < 0) {
        st = h5c__fail_hdf5(-1, "H5Awrite failed for '%s' on '%s'",
                            name, obj_path);
    }

done:
    if (aid >= 0) {
        H5Aclose(aid);
    }
    if (sid >= 0) {
        H5Sclose(sid);
    }
    H5Tclose(tid);
    free(buf);
    H5Oclose(oid);
    return st;
}

h5c_status_t h5c_write_attr_str(h5c_file_t *file, const char *obj_path,
                                const char *name, const char *value)
{
    return h5c__record(file, write_str_impl(file, obj_path, name, value));
}

static h5c_status_t read_str_impl(h5c_file_t *file, const char *obj_path,
                                  const char *name, char **out)
{
    h5c_status_t st;
    hid_t        oid = H5I_INVALID_HID, aid = H5I_INVALID_HID;

    if (out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_read_attr_str: out is NULL");
    }
    *out = NULL;

    if ((st = check_args(file, obj_path, name)) != H5C_OK) {
        return st;
    }
    if ((st = open_object(file, obj_path, &oid)) != H5C_OK) {
        return st;
    }
    if ((st = open_attr(oid, obj_path, name, &aid)) != H5C_OK) {
        H5Oclose(oid);
        return st;
    }

    st = h5c__str_read_id(aid, 1, name, out);

    H5Aclose(aid);
    H5Oclose(oid);
    return st;
}

h5c_status_t h5c_read_attr_str(h5c_file_t *file, const char *obj_path,
                               const char *name, char **out)
{
    return h5c__record(file, read_str_impl(file, obj_path, name, out));
}

/* ------------------------------------------------------------------ */
/* scalar numeric attributes                                           */
/* ------------------------------------------------------------------ */

static h5c_status_t write_scalar_impl(h5c_file_t *file, const char *obj_path,
                                      const char *name, const void *value,
                                      h5c_type_t type)
{
    h5c_status_t st;
    hid_t        oid = H5I_INVALID_HID, sid = H5I_INVALID_HID;
    hid_t        aid = H5I_INVALID_HID, ftype, mtype;

    if ((st = check_args(file, obj_path, name)) != H5C_OK) {
        return st;
    }
    if (value == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "value is NULL for attribute '%s'", name);
    }
    if (file->readonly) {
        return h5c__fail(H5C_ERR_STATE,
                         "file is open read-only, cannot write '%s'", name);
    }

    ftype = h5c__file_type(type);
    mtype = h5c__mem_type(type);
    if (ftype == H5I_INVALID_HID || mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping "
                         "(use h5c_write_attr_str for strings)", (int)type);
    }
    if ((st = open_object(file, obj_path, &oid)) != H5C_OK) {
        return st;
    }

    sid = H5Screate(H5S_SCALAR);
    if (sid < 0) {
        H5Oclose(oid);
        return h5c__fail_hdf5((long)sid,
                              "cannot create a scalar space for '%s'", name);
    }
    if ((st = drop_existing(oid, obj_path, name)) != H5C_OK) {
        goto done;
    }

    aid = H5Acreate2(oid, name, ftype, sid, H5P_DEFAULT, H5P_DEFAULT);
    if (aid < 0) {
        st = h5c__fail_hdf5((long)aid, "cannot create attribute '%s' on '%s'",
                            name, obj_path);
        goto done;
    }
    if (H5Awrite(aid, mtype, value) < 0) {
        st = h5c__fail_hdf5(-1, "H5Awrite failed for '%s' on '%s'",
                            name, obj_path);
    }

done:
    if (aid >= 0) {
        H5Aclose(aid);
    }
    H5Sclose(sid);
    H5Oclose(oid);
    return st;
}

h5c_status_t h5c_write_attr_scalar(h5c_file_t *file, const char *obj_path,
                                   const char *name, const void *value,
                                   h5c_type_t type)
{
    return h5c__record(file,
                       write_scalar_impl(file, obj_path, name, value, type));
}

static h5c_status_t read_scalar_impl(h5c_file_t *file, const char *obj_path,
                                     const char *name, void *value,
                                     h5c_type_t type)
{
    h5c_status_t st;
    hid_t        oid = H5I_INVALID_HID, aid = H5I_INVALID_HID;
    hid_t        sid = H5I_INVALID_HID, mtype;

    if ((st = check_args(file, obj_path, name)) != H5C_OK) {
        return st;
    }
    if (value == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "value is NULL for attribute '%s'", name);
    }

    mtype = h5c__mem_type_read(type);
    if (mtype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping "
                         "(use h5c_read_attr_str for strings)", (int)type);
    }
    if ((st = open_object(file, obj_path, &oid)) != H5C_OK) {
        return st;
    }
    if ((st = open_attr(oid, obj_path, name, &aid)) != H5C_OK) {
        H5Oclose(oid);
        return st;
    }

    sid = H5Aget_space(aid);
    if (sid < 0) {
        st = h5c__fail_hdf5((long)sid, "cannot query the space of '%s'", name);
        goto done;
    }
    if (H5Sget_simple_extent_npoints(sid) != 1) {
        st = h5c__fail(H5C_ERR_SHAPE_MISMATCH,
                       "attribute '%s' on '%s' is not a scalar",
                       name, obj_path);
        goto done;
    }
    /* HDF5 converts between the stored type and `mtype` where it can. */
    if (H5Aread(aid, mtype, value) < 0) {
        st = h5c__fail_hdf5(-1, "H5Aread failed for '%s' on '%s'",
                            name, obj_path);
    }

done:
    if (sid >= 0) {
        H5Sclose(sid);
    }
    H5Aclose(aid);
    H5Oclose(oid);
    return st;
}

h5c_status_t h5c_read_attr_scalar(h5c_file_t *file, const char *obj_path,
                                  const char *name, void *value,
                                  h5c_type_t type)
{
    return h5c__record(file,
                       read_scalar_impl(file, obj_path, name, value, type));
}

/* ------------------------------------------------------------------ */
/* existence                                                           */
/* ------------------------------------------------------------------ */

/* Like h5c_exists(), this is a query: it never touches the sticky status. */
int h5c_attr_exists(h5c_file_t *file, const char *obj_path, const char *name)
{
    hid_t  oid;
    htri_t present;

    if (file == NULL || file->fid < 0 ||
        obj_path == NULL || obj_path[0] == '\0' ||
        name == NULL || name[0] == '\0') {
        return 0;
    }
    if (h5c__ensure_init() != H5C_OK) {
        return 0;
    }
    oid = H5Oopen(file->fid, obj_path, H5P_DEFAULT);
    if (oid < 0) {
        return 0;
    }
    present = H5Aexists(oid, name);
    H5Oclose(oid);
    return (present > 0) ? 1 : 0;
}
