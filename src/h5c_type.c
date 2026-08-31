#include "h5c_internal.h"

/*
 * Cached ids that h5c owns for the whole process:
 *   - the boolean enum, used as BOTH the file and the memory type so that
 *     writing an int8 buffer is a plain copy with no conversion path;
 *   - a link creation property list that creates intermediate groups.
 */
static hid_t g_bool_type = H5I_INVALID_HID;
static hid_t g_lcpl      = H5I_INVALID_HID;

static hid_t bool_type(void)
{
    if (g_bool_type == H5I_INVALID_HID) {
        const int8_t f = 0, t = 1;
        hid_t tid = H5Tenum_create(H5T_STD_I8LE);
        if (tid < 0) {
            h5c__fail_hdf5((long)tid, "H5Tenum_create failed for bool");
            return H5I_INVALID_HID;
        }
        if (H5Tenum_insert(tid, "FALSE", &f) < 0 ||
            H5Tenum_insert(tid, "TRUE",  &t) < 0) {
            H5Tclose(tid);
            h5c__fail_hdf5(-1, "H5Tenum_insert failed for bool");
            return H5I_INVALID_HID;
        }
        g_bool_type = tid;
    }
    return g_bool_type;
}

hid_t h5c__file_type(h5c_type_t type)
{
    switch (type) {
    /* Explicit little-endian so files are reproducible across platforms. */
    case H5C_F32:  return H5T_IEEE_F32LE;
    case H5C_F64:  return H5T_IEEE_F64LE;
    case H5C_I32:  return H5T_STD_I32LE;
    case H5C_I64:  return H5T_STD_I64LE;
    case H5C_BOOL: return bool_type();
    default:       return H5I_INVALID_HID;
    }
}

hid_t h5c__mem_type(h5c_type_t type)
{
    switch (type) {
    case H5C_F32:  return H5T_NATIVE_FLOAT;
    case H5C_F64:  return H5T_NATIVE_DOUBLE;
    case H5C_I32:  return H5T_NATIVE_INT32;
    case H5C_I64:  return H5T_NATIVE_INT64;
    case H5C_BOOL: return bool_type();
    default:       return H5I_INVALID_HID;
    }
}

size_t h5c_type_size(h5c_type_t type)
{
    switch (type) {
    case H5C_F32:  return sizeof(float);
    case H5C_F64:  return sizeof(double);
    case H5C_I32:  return sizeof(int32_t);
    case H5C_I64:  return sizeof(int64_t);
    case H5C_BOOL: return sizeof(h5c_bool_t);
    default:       return 0;
    }
}

h5c_type_t h5c__type_from_hid(hid_t tid)
{
    H5T_class_t cls = H5Tget_class(tid);
    size_t      sz  = H5Tget_size(tid);

    switch (cls) {
    case H5T_FLOAT:
        if (sz == sizeof(float))  return H5C_F32;
        if (sz == sizeof(double)) return H5C_F64;
        return H5C_TYPE_UNKNOWN;
    case H5T_INTEGER:
        if (sz == 8) return H5C_I64;
        if (sz == 4) return H5C_I32;
        /* One-byte integers are how a plain (non-enum) bool would be stored. */
        if (sz == 1) return H5C_BOOL;
        return H5C_TYPE_UNKNOWN;
    case H5T_ENUM:
        return H5C_BOOL;
    case H5T_STRING:
        return H5C_STRING;
    default:
        return H5C_TYPE_UNKNOWN;
    }
}

hid_t h5c__lcpl(void)
{
    if (g_lcpl == H5I_INVALID_HID) {
        hid_t p = H5Pcreate(H5P_LINK_CREATE);
        if (p < 0) {
            h5c__fail_hdf5((long)p, "H5Pcreate(H5P_LINK_CREATE) failed");
            return H5P_DEFAULT;
        }
        if (H5Pset_create_intermediate_group(p, 1) < 0) {
            H5Pclose(p);
            h5c__fail_hdf5(-1, "H5Pset_create_intermediate_group failed");
            return H5P_DEFAULT;
        }
        g_lcpl = p;
    }
    return g_lcpl;
}

void h5c__type_cleanup(void)
{
    if (g_bool_type != H5I_INVALID_HID) {
        H5Tclose(g_bool_type);
        g_bool_type = H5I_INVALID_HID;
    }
    if (g_lcpl != H5I_INVALID_HID) {
        H5Pclose(g_lcpl);
        g_lcpl = H5I_INVALID_HID;
    }
}
