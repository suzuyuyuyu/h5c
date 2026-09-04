/*
 * Cross-language format regression tests.
 *
 * This file is the executable counterpart of docs/FORMAT.md. It pins down the
 * on-disk layout that h5c shares with h5fortran, in both directions:
 *
 *   REFERENCE WRITER -> h5c reader
 *       Builds a file with the BARE HDF5 C API only (H5Fcreate, H5Screate_simple,
 *       H5Dcreate2, ...). It never calls h5c. h5c then reads it and must see the
 *       declared shapes, types and values.
 *
 *   h5c writer -> REFERENCE READER
 *       h5c writes a file, the file is reopened with H5Fopen, and the raw
 *       dataspace extents, datatype properties and bytes are inspected
 *       directly. A regression is then reported as "dims are reversed" rather
 *       than merely "values differ".
 *
 * No reference .h5 is committed: binary blobs have unreadable diffs. The
 * expected values live here, in source, where they are diffable.
 *
 * Every shape and every value is ASYMMETRIC on purpose. A transposed or
 * dimension-reversing implementation must not be able to reproduce them by
 * accident, so there is no 2x2 array and no symmetric data anywhere below.
 */
#include "h5c_test.h"

#include <stdlib.h>
#include <string.h>

H5C_TEST_MAIN_STATE;

#define REF_PATH "test_crosslang_ref.h5"   /* written by the bare HDF5 API */
#define H5C_PATH "test_crosslang_h5c.h5"   /* written by h5c               */

#define H5FORTRAN_ARTIFACT \
    "/home/b/b39007/workspace/dev-hdf5/h5fortran/build-integration/test/test-serial.h5"

/* ------------------------------------------------------------------ */
/* the shared expectation                                              */
/* ------------------------------------------------------------------ */

/*
 * The canonical case from docs/FORMAT.md: a Fortran a(2,3) holding 1..6 is
 * stored as dataspace {3, 2} with the flat byte sequence 1 2 3 4 5 6.
 */
static const double k_2d[6]        = { 1, 2, 3, 4, 5, 6 };
static const size_t k_2d_dims[2]   = { 3, 2 };
static const hsize_t k_2d_hdims[2] = { 3, 2 };

/*
 * Rank 3, no two extents equal, and each value encodes its own HDF5
 * coordinate as 100*i0 + 10*i1 + i2. Any permutation of the dimensions is
 * therefore detectable from a single element.
 */
#define K_3D_D0 4u
#define K_3D_D1 3u
#define K_3D_D2 2u
#define K_3D_N  (K_3D_D0 * K_3D_D1 * K_3D_D2)

static const size_t  k_3d_dims[3]  = { K_3D_D0, K_3D_D1, K_3D_D2 };
static const hsize_t k_3d_hdims[3] = { K_3D_D0, K_3D_D1, K_3D_D2 };

static void fill_3d(double *buf)
{
    unsigned i0, i1, i2;

    for (i0 = 0; i0 < K_3D_D0; i0++) {
        for (i1 = 0; i1 < K_3D_D1; i1++) {
            for (i2 = 0; i2 < K_3D_D2; i2++) {
                /* row-major: dims[rank-1] varies fastest */
                buf[(i0 * K_3D_D1 + i1) * K_3D_D2 + i2] =
                    (double)(100u * i0 + 10u * i1 + i2);
            }
        }
    }
}

/* Only the raw-byte comparisons depend on this; the value checks do not. */
static int host_is_le(void)
{
    const uint32_t one = 1u;
    unsigned char  b[sizeof one];

    memcpy(b, &one, sizeof one);
    return b[0] == 1u;
}

/* ------------------------------------------------------------------ */
/* bare-HDF5 helpers (no h5c call in this section)                     */
/* ------------------------------------------------------------------ */

static hid_t g_ref_lcpl = H5I_INVALID_HID;

static hid_t ref_lcpl(void)
{
    if (g_ref_lcpl == H5I_INVALID_HID) {
        hid_t p = H5Pcreate(H5P_LINK_CREATE);

        if (p < 0 || H5Pset_create_intermediate_group(p, 1) < 0) {
            H5C_FAILF("cannot build a link creation property list");
            return H5P_DEFAULT;
        }
        g_ref_lcpl = p;
    }
    return g_ref_lcpl;
}

/*
 * Creates `path` with the given file type and writes `buf` through `mtype`.
 * rank == 0 makes a scalar dataset. `buf` may be NULL for an empty extent.
 */
static void ref_write(hid_t fid, const char *path, hid_t ftype, hid_t mtype,
                      int rank, const hsize_t *dims, const void *buf)
{
    hid_t     sid;
    hid_t     did;
    hssize_t  npoints;

    sid = (rank == 0) ? H5Screate(H5S_SCALAR)
                      : H5Screate_simple(rank, dims, NULL);
    if (sid < 0) {
        H5C_FAILF("H5Screate_simple failed for '%s'", path);
        return;
    }
    did = H5Dcreate2(fid, path, ftype, sid, ref_lcpl(), H5P_DEFAULT,
                     H5P_DEFAULT);
    if (did < 0) {
        H5C_FAILF("H5Dcreate2 failed for '%s'", path);
        H5Sclose(sid);
        return;
    }
    npoints = H5Sget_simple_extent_npoints(sid);
    if (npoints > 0 && buf != NULL) {
        if (H5Dwrite(did, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) {
            H5C_FAILF("H5Dwrite failed for '%s'", path);
        }
    }
    H5Dclose(did);
    H5Sclose(sid);
}

static void ref_write_attr(hid_t fid, const char *obj_path, const char *name,
                           hid_t ftype, hid_t mtype, hsize_t count,
                           const void *buf)
{
    hid_t oid, sid, aid;

    oid = H5Oopen(fid, obj_path, H5P_DEFAULT);
    if (oid < 0) {
        H5C_FAILF("H5Oopen failed for '%s'", obj_path);
        return;
    }
    sid = H5Screate_simple(1, &count, NULL);
    if (sid < 0) {
        H5C_FAILF("H5Screate_simple failed for attribute '%s'", name);
        H5Oclose(oid);
        return;
    }
    aid = H5Acreate2(oid, name, ftype, sid, H5P_DEFAULT, H5P_DEFAULT);
    if (aid < 0) {
        H5C_FAILF("H5Acreate2 failed for attribute '%s'", name);
    } else {
        if (H5Awrite(aid, mtype, buf) < 0) {
            H5C_FAILF("H5Awrite failed for attribute '%s'", name);
        }
        H5Aclose(aid);
    }
    H5Sclose(sid);
    H5Oclose(oid);
}

static void ref_check_attr_dims(hid_t fid, const char *obj_path,
                                const char *name, hsize_t want)
{
    hid_t   oid, aid, sid;
    hsize_t got[1] = { 0 };
    int     rank;

    oid = H5Oopen(fid, obj_path, H5P_DEFAULT);
    aid = (oid >= 0) ? H5Aopen(oid, name, H5P_DEFAULT) : H5I_INVALID_HID;
    sid = (aid >= 0) ? H5Aget_space(aid) : H5I_INVALID_HID;
    if (oid < 0 || aid < 0 || sid < 0) {
        H5C_FAILF("cannot inspect attribute '%s' on '%s'", name, obj_path);
    } else {
        rank = H5Sget_simple_extent_ndims(sid);
        H5C_ASSERT(rank == 1, "attribute '%s' rank %d, want 1", name, rank);
        if (rank == 1 && H5Sget_simple_extent_dims(sid, got, NULL) >= 0) {
            H5C_ASSERT_EQ_SIZE(got[0], want, name);
        } else if (rank == 1) {
            H5C_FAILF("H5Sget_simple_extent_dims failed for '%s'", name);
        }
    }
    if (sid >= 0) {
        H5Sclose(sid);
    }
    if (aid >= 0) {
        H5Aclose(aid);
    }
    if (oid >= 0) {
        H5Oclose(oid);
    }
}

/* The boolean enum as docs/FORMAT.md declares it. The caller closes it. */
static hid_t ref_bool_type(void)
{
    const int8_t f = 0, t = 1;
    hid_t        tid = H5Tenum_create(H5T_STD_I8LE);

    if (tid < 0 || H5Tenum_insert(tid, "FALSE", &f) < 0 ||
        H5Tenum_insert(tid, "TRUE", &t) < 0) {
        H5C_FAILF("cannot build the reference boolean enum");
        return H5I_INVALID_HID;
    }
    return tid;
}

/* Fixed-length H5T_C_S1, SPACEPAD, ASCII. The caller closes it. */
static hid_t ref_string_type(size_t len)
{
    hid_t tid = H5Tcopy(H5T_C_S1);

    if (tid < 0 || H5Tset_size(tid, len) < 0 ||
        H5Tset_strpad(tid, H5T_STR_SPACEPAD) < 0 ||
        H5Tset_cset(tid, H5T_CSET_ASCII) < 0) {
        H5C_FAILF("cannot build the reference string type");
        return H5I_INVALID_HID;
    }
    return tid;
}

/*
 * Asserts the dataspace extents of `path` element by element. This is the
 * check that names a reversed-dimension regression outright, because it
 * compares H5Sget_simple_extent_dims() against the expected order instead of
 * comparing values that a transposed square case would satisfy anyway.
 */
static void ref_check_dims(hid_t fid, const char *path, int want_rank,
                           const hsize_t *want_dims)
{
    hid_t   did, sid;
    hsize_t got[H5C_MAX_RANK];
    int     rank, i;

    did = H5Dopen2(fid, path, H5P_DEFAULT);
    if (did < 0) {
        H5C_FAILF("H5Dopen2 failed for '%s'", path);
        return;
    }
    sid = H5Dget_space(did);
    if (sid < 0) {
        H5C_FAILF("H5Dget_space failed for '%s'", path);
        H5Dclose(did);
        return;
    }
    rank = H5Sget_simple_extent_ndims(sid);
    if (rank != want_rank) {
        H5C_FAILF("'%s': rank %d, want %d", path, rank, want_rank);
    } else if (rank > 0) {
        if (H5Sget_simple_extent_dims(sid, got, NULL) < 0) {
            H5C_FAILF("H5Sget_simple_extent_dims failed for '%s'", path);
        } else {
            for (i = 0; i < rank; i++) {
                if (got[i] != want_dims[i]) {
                    H5C_FAILF("'%s': dims[%d] is %lu, want %lu"
                              " (a reversed-dimension regression looks"
                              " exactly like this)",
                              path, i, (unsigned long)got[i],
                              (unsigned long)want_dims[i]);
                }
            }
        }
    }
    H5Sclose(sid);
    H5Dclose(did);
}

/* Asserts the stored datatype's class, size, signedness and byte order. */
static void ref_check_numeric_type(hid_t fid, const char *path,
                                   H5T_class_t want_class, size_t want_size,
                                   int want_signed)
{
    hid_t       did, tid;
    H5T_class_t cls;
    size_t      sz;
    H5T_order_t order;

    did = H5Dopen2(fid, path, H5P_DEFAULT);
    if (did < 0) {
        H5C_FAILF("H5Dopen2 failed for '%s'", path);
        return;
    }
    tid = H5Dget_type(did);
    if (tid < 0) {
        H5C_FAILF("H5Dget_type failed for '%s'", path);
        H5Dclose(did);
        return;
    }
    cls   = H5Tget_class(tid);
    sz    = H5Tget_size(tid);
    order = H5Tget_order(tid);

    H5C_ASSERT(cls == want_class, "'%s': datatype class %d, want %d", path,
               (int)cls, (int)want_class);
    H5C_ASSERT_EQ_SIZE(sz, want_size, path);
    /* A one-byte type has no meaningful byte order. */
    if (sz > 1) {
        H5C_ASSERT(order == H5T_ORDER_LE,
                   "'%s': byte order %d, want little-endian (%d)", path,
                   (int)order, (int)H5T_ORDER_LE);
    }
    if (want_signed && cls == H5T_INTEGER) {
        H5C_ASSERT(H5Tget_sign(tid) == H5T_SGN_2,
                   "'%s': integer is not two's-complement signed", path);
    }
    H5Tclose(tid);
    H5Dclose(did);
}

/*
 * Reads the raw file bytes of `path` by asking for the stored type itself, so
 * no conversion can hide a layout error, and compares them with `want`.
 */
static void ref_check_bytes(hid_t fid, const char *path,
                            const void *want, size_t nbytes)
{
    hid_t          did, tid;
    unsigned char *got;

    if (!host_is_le()) {
        return; /* the value checks still cover this case */
    }
    did = H5Dopen2(fid, path, H5P_DEFAULT);
    if (did < 0) {
        H5C_FAILF("H5Dopen2 failed for '%s'", path);
        return;
    }
    tid = H5Dget_type(did);
    got = (unsigned char *)calloc(nbytes > 0 ? nbytes : 1, 1);
    if (tid < 0 || got == NULL) {
        H5C_FAILF("cannot stage a raw read of '%s'", path);
        free(got);
        if (tid >= 0) {
            H5Tclose(tid);
        }
        H5Dclose(did);
        return;
    }
    if (H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, got) < 0) {
        H5C_FAILF("raw H5Dread failed for '%s'", path);
    } else if (memcmp(got, want, nbytes) != 0) {
        size_t i;

        H5C_FAILF("'%s': raw bytes differ from the expected file image", path);
        for (i = 0; i < nbytes; i++) {
            if (got[i] != ((const unsigned char *)want)[i]) {
                fprintf(stderr, "       first difference at byte %lu:"
                                " got 0x%02x, want 0x%02x\n",
                        (unsigned long)i, got[i],
                        ((const unsigned char *)want)[i]);
                break;
            }
        }
    }
    free(got);
    H5Tclose(tid);
    H5Dclose(did);
}

/* ------------------------------------------------------------------ */
/* phase A: the reference writer                                       */
/* ------------------------------------------------------------------ */

static void build_reference_file(void)
{
    hid_t   fid;
    hid_t   bt, st;
    double  d3[K_3D_N];
    /* h5fortran writes `logical` as H5T_STD_I32LE holding 1/0. */
    const int32_t lgc_i32[6] = { 1, 0, 0, 1, 1, 0 };
    const int8_t  flags[6]   = { 1, 0, 0, 1, 1, 0 };
    const int32_t i32[6]     = { 10, -20, 30, -40, 50, -60 };
    const int64_t i64[6]     = { 1, -2, 3000000000LL, -4000000000LL, 5, -6 };
    const int16_t i16[6]     = { 1, -2, 3, -4, 5, -6 };
    const int8_t  i8[6]      = { 1, -2, 3, -4, 5, -6 };
    const float   f32[6]     = { 1.5f, -2.5f, 3.5f, -4.5f, 5.5f, -6.5f };
    /* Fixed length 16, space padded: what h5c must trim back to "padded". */
    const char    padded[16] = { 'p', 'a', 'd', 'd', 'e', 'd',
                                 ' ', ' ', ' ', ' ', ' ', ' ',
                                 ' ', ' ', ' ', ' ' };
    const hsize_t zero_dims[2] = { 0, 2 };
    const int32_t attr_i32[4] = { 71, -8, 3005, 42 };
    const int32_t attr_bool_i32[4] = { 1, 1, 0, 1 };

    fill_3d(d3);

    fid = H5Fcreate(REF_PATH, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fid < 0) {
        H5C_FAILF("H5Fcreate failed for " REF_PATH);
        return;
    }

    ref_write(fid, "/rank/two", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 2,
              k_2d_hdims, k_2d);
    ref_write_attr(fid, "/rank/two", "reference_array", H5T_STD_I32LE,
                   H5T_NATIVE_INT32, 4, attr_i32);
    ref_check_attr_dims(fid, "/rank/two", "reference_array", 4);
    ref_write_attr(fid, "/rank/two", "reference_bool_array", H5T_STD_I32LE,
                   H5T_NATIVE_INT32, 4, attr_bool_i32);
    ref_write(fid, "/rank/three", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 3,
              k_3d_hdims, d3);
    ref_write(fid, "/scalar/r64", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 0, NULL,
              &k_2d[3]);

    ref_write(fid, "/types/f32", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, 2,
              k_2d_hdims, f32);
    ref_write(fid, "/types/i8", H5T_STD_I8LE, H5T_NATIVE_INT8, 2,
              k_2d_hdims, i8);
    ref_write(fid, "/types/i16", H5T_STD_I16LE, H5T_NATIVE_INT16, 2,
              k_2d_hdims, i16);
    ref_write(fid, "/types/i32", H5T_STD_I32LE, H5T_NATIVE_INT32, 2,
              k_2d_hdims, i32);
    ref_write(fid, "/types/i64", H5T_STD_I64LE, H5T_NATIVE_INT64, 2,
              k_2d_hdims, i64);

    bt = ref_bool_type();
    if (bt >= 0) {
        ref_write(fid, "/types/bool_enum", bt, bt, 2, k_2d_hdims, flags);
        H5Tclose(bt);
    }
    /* The h5fortran representation, which h5c must also accept. */
    ref_write(fid, "/types/bool_i32", H5T_STD_I32LE, H5T_NATIVE_INT32, 2,
              k_2d_hdims, lgc_i32);

    st = ref_string_type(sizeof padded);
    if (st >= 0) {
        ref_write(fid, "/types/text", st, st, 0, NULL, padded);
        H5Tclose(st);
    }

    /* A zero extent is a real object with a declared shape. */
    ref_write(fid, "/empty/two", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 2,
              zero_dims, NULL);

    if (H5Fclose(fid) < 0) {
        H5C_FAILF("H5Fclose failed for " REF_PATH);
    }
}

/* ------------------------------------------------------------------ */
/* phase B: h5c reads what the reference writer produced               */
/* ------------------------------------------------------------------ */

/*
 * KNOWN GAP, measured here rather than assumed.
 *
 * h5fortran writes `logical` as an H5T_STD_I32LE 1/0 dataset, and h5c must be
 * able to read it as H5C_BOOL. That is not automatic: HDF5 converts an enum to
 * an integer but has NO integer -> enum path, so reading through the boolean
 * enum fails on an integer dataset. h5c therefore reads booleans through
 * H5T_NATIVE_INT8 (h5c__mem_type_read) and only writes through the enum.
 *
 * This test guards that asymmetry, which an innocent-looking "use one memory
 * type everywhere" cleanup would undo.
 */
static void check_int32_logical(h5c_file_t *f, const char *path,
                                const h5c_bool_t *want, size_t n,
                                const size_t *dims)
{
    h5c_bool_t   *gb = (h5c_bool_t *)calloc(n, sizeof *gb);
    int8_t       *g8 = (int8_t *)calloc(n, sizeof *g8);
    h5c_status_t  st;
    size_t        i;

    if (gb == NULL || g8 == NULL) {
        H5C_FAILF("out of memory");
        free(gb);
        free(g8);
        return;
    }

    /* The interoperability that matters: h5fortran's logical read as bool. */
    st = h5c_read_bool(f, path, gb, 2, dims);
    H5C_ASSERT(st == H5C_OK, "h5c_read_bool('%s') -> %s", path,
               h5c_status_string(st));
    if (st == H5C_OK) {
        for (i = 0; i < n; i++) {
            H5C_ASSERT(gb[i] == want[i], "'%s' as bool[%lu]: got %d want %d",
                       path, (unsigned long)i, (int)gb[i], (int)want[i]);
        }
    }

    /* This path does work, and it is what makes the data recoverable. */
    H5C_CHECK(h5c_read_i8(f, path, g8, 2, dims));
    for (i = 0; i < n; i++) {
        H5C_ASSERT(g8[i] == (int8_t)want[i], "'%s' as i8[%lu]: got %d want %d",
                   path, (unsigned long)i, (int)g8[i], (int)want[i]);
    }

    free(gb);
    free(g8);
}

static void read_reference_shapes(h5c_file_t *f)
{
    h5c_dataset_info_t info;
    double             got2[6];
    double             got3[K_3D_N];
    unsigned           i0, i1, i2;
    int                i;

    /* rank 2: dims must be reported as {3, 2}, in that order. */
    H5C_CHECK(h5c_dataset_info(f, "/rank/two", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 2, "rank of the reference /rank/two");
    H5C_ASSERT_EQ_SIZE(info.dims[0], 3, "dims[0] of the reference /rank/two");
    H5C_ASSERT_EQ_SIZE(info.dims[1], 2, "dims[1] of the reference /rank/two");
    H5C_ASSERT_EQ_SIZE(info.count, 6, "count of the reference /rank/two");
    H5C_ASSERT(info.type == H5C_F64, "type of the reference /rank/two: %d",
               (int)info.type);

    memset(got2, 0, sizeof got2);
    H5C_CHECK(h5c_read_f64(f, "/rank/two", got2, 2, k_2d_dims));
    for (i = 0; i < 6; i++) {
        H5C_ASSERT(got2[i] == k_2d[i], "reference 2d flat[%d]: got %g want %g",
                   i, got2[i], k_2d[i]);
    }

    /*
     * Swapping the dimensions must be REJECTED, never silently transposed.
     * With a {3,2} dataset a {2,3} request has the same element count, so
     * nothing but an explicit extent comparison can catch it.
     */
    {
        const size_t swapped[2] = { 2, 3 };

        H5C_CHECK_FAILS(h5c_read_f64(f, "/rank/two", got2, 2, swapped),
                        H5C_ERR_SHAPE_MISMATCH);
    }

    /* rank 3: every value carries its coordinate, so any permutation shows. */
    H5C_CHECK(h5c_dataset_info(f, "/rank/three", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 3, "rank of the reference /rank/three");
    H5C_ASSERT_EQ_SIZE(info.dims[0], K_3D_D0, "dims[0] of /rank/three");
    H5C_ASSERT_EQ_SIZE(info.dims[1], K_3D_D1, "dims[1] of /rank/three");
    H5C_ASSERT_EQ_SIZE(info.dims[2], K_3D_D2, "dims[2] of /rank/three");

    memset(got3, 0, sizeof got3);
    H5C_CHECK(h5c_read_f64(f, "/rank/three", got3, 3, k_3d_dims));
    for (i0 = 0; i0 < K_3D_D0; i0++) {
        for (i1 = 0; i1 < K_3D_D1; i1++) {
            for (i2 = 0; i2 < K_3D_D2; i2++) {
                size_t flat = (i0 * K_3D_D1 + i1) * K_3D_D2 + i2;
                double want = (double)(100u * i0 + 10u * i1 + i2);

                H5C_ASSERT(got3[flat] == want,
                           "reference 3d (%u,%u,%u) at flat %lu:"
                           " got %g want %g",
                           i0, i1, i2, (unsigned long)flat, got3[flat], want);
            }
        }
    }

    /* Both non-trivial permutations of {4,3,2} must be refused. */
    {
        const size_t perm_a[3] = { K_3D_D2, K_3D_D1, K_3D_D0 }; /* 2,3,4 */
        const size_t perm_b[3] = { K_3D_D1, K_3D_D0, K_3D_D2 }; /* 3,4,2 */

        H5C_CHECK_FAILS(h5c_read_f64(f, "/rank/three", got3, 3, perm_a),
                        H5C_ERR_SHAPE_MISMATCH);
        H5C_CHECK_FAILS(h5c_read_f64(f, "/rank/three", got3, 3, perm_b),
                        H5C_ERR_SHAPE_MISMATCH);
    }

    /* Scalars have rank 0 and count 1. */
    H5C_CHECK(h5c_dataset_info(f, "/scalar/r64", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 0, "rank of the reference /scalar/r64");
    H5C_ASSERT_EQ_SIZE(info.count, 1, "count of the reference /scalar/r64");

    /* The provoked shape failures are on the sticky record by design. */
    h5c_file_clear_status(f);
}

static void read_reference_types(h5c_file_t *f)
{
    h5c_dataset_info_t info;
    const int32_t      want_i32[6] = { 10, -20, 30, -40, 50, -60 };
    const int64_t      want_i64[6] = { 1, -2, 3000000000LL, -4000000000LL,
                                       5, -6 };
    const int16_t      want_i16[6] = { 1, -2, 3, -4, 5, -6 };
    const int8_t       want_i8[6]  = { 1, -2, 3, -4, 5, -6 };
    const float        want_f32[6] = { 1.5f, -2.5f, 3.5f, -4.5f, 5.5f, -6.5f };
    const h5c_bool_t   want_b[6]   = { H5C_TRUE, H5C_FALSE, H5C_FALSE,
                                       H5C_TRUE, H5C_TRUE, H5C_FALSE };
    int32_t     g32[6];
    int64_t     g64[6];
    int16_t     g16[6];
    int8_t      g8[6];
    float       gf32[6];
    h5c_bool_t  gb[6];
    char       *text = NULL;
    int         i;

    H5C_CHECK(h5c_read_f32(f, "/types/f32", gf32, 2, k_2d_dims));
    H5C_CHECK(h5c_read_i8(f, "/types/i8", g8, 2, k_2d_dims));
    H5C_CHECK(h5c_read_i16(f, "/types/i16", g16, 2, k_2d_dims));
    H5C_CHECK(h5c_read_i32(f, "/types/i32", g32, 2, k_2d_dims));
    H5C_CHECK(h5c_read_i64(f, "/types/i64", g64, 2, k_2d_dims));
    for (i = 0; i < 6; i++) {
        H5C_ASSERT(gf32[i] == want_f32[i], "reference f32[%d]: got %g", i,
                   (double)gf32[i]);
        H5C_ASSERT(g8[i] == want_i8[i], "reference i8[%d]: got %d", i,
                   (int)g8[i]);
        H5C_ASSERT(g16[i] == want_i16[i], "reference i16[%d]: got %d", i,
                   (int)g16[i]);
        H5C_ASSERT(g32[i] == want_i32[i], "reference i32[%d]: got %d", i,
                   (int)g32[i]);
        H5C_ASSERT(g64[i] == want_i64[i], "reference i64[%d] mismatch", i);
    }

    /* The type that h5c reports back for each stored type. */
    H5C_CHECK(h5c_dataset_info(f, "/types/f32", &info));
    H5C_ASSERT(info.type == H5C_F32, "/types/f32 reports %d", (int)info.type);
    H5C_CHECK(h5c_dataset_info(f, "/types/i8", &info));
    H5C_ASSERT(info.type == H5C_I8, "/types/i8 reports %d", (int)info.type);
    H5C_CHECK(h5c_dataset_info(f, "/types/i16", &info));
    H5C_ASSERT(info.type == H5C_I16, "/types/i16 reports %d", (int)info.type);
    H5C_CHECK(h5c_dataset_info(f, "/types/i32", &info));
    H5C_ASSERT(info.type == H5C_I32, "/types/i32 reports %d", (int)info.type);
    H5C_CHECK(h5c_dataset_info(f, "/types/i64", &info));
    H5C_ASSERT(info.type == H5C_I64, "/types/i64 reports %d", (int)info.type);
    H5C_CHECK(h5c_dataset_info(f, "/types/bool_enum", &info));
    H5C_ASSERT(info.type == H5C_BOOL, "/types/bool_enum reports %d",
               (int)info.type);

    /* The enum form: no conversion, a plain copy. */
    memset(gb, 0x7f, sizeof gb);
    H5C_CHECK(h5c_read_bool(f, "/types/bool_enum", gb, 2, k_2d_dims));
    for (i = 0; i < 6; i++) {
        H5C_ASSERT(gb[i] == want_b[i], "reference bool_enum[%d]: got %d", i,
                   (int)gb[i]);
    }

    /* The int32 0/1 form, which is what h5fortran writes for `logical`. */
    check_int32_logical(f, "/types/bool_i32", want_b, 6, k_2d_dims);

    /* SPACEPAD means the trailing spaces are padding, not content. */
    H5C_CHECK(h5c_read_string(f, "/types/text", &text));
    if (text != NULL) {
        H5C_ASSERT(strcmp(text, "padded") == 0,
                   "space-padded string: got '%s', want 'padded'", text);
        h5c_free_string(text);
    }
}

static void read_reference_empty(h5c_file_t *f)
{
    h5c_dataset_info_t info;
    const size_t       dims[2] = { 0, 2 };
    double            *alloc   = NULL;

    H5C_ASSERT(h5c_exists(f, "/empty/two"),
               "an empty dataset is still a real object");
    H5C_CHECK(h5c_dataset_info(f, "/empty/two", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 2, "rank of the reference /empty/two");
    H5C_ASSERT_EQ_SIZE(info.dims[0], 0, "dims[0] of the reference /empty/two");
    H5C_ASSERT_EQ_SIZE(info.dims[1], 2, "dims[1] of the reference /empty/two");
    H5C_ASSERT_EQ_SIZE(info.count, 0, "count of the reference /empty/two");

    /* Zero elements: the read succeeds and touches nothing. */
    H5C_CHECK(h5c_read_f64(f, "/empty/two", NULL, 2, dims));

    H5C_CHECK(h5c_read_alloc(f, "/empty/two", H5C_F64, (void **)&alloc, &info));
    H5C_ASSERT(alloc != NULL, "read_alloc stays non-NULL for an empty dataset");
    h5c_free(alloc);
}

static void read_reference_attr(h5c_file_t *f)
{
    const int32_t want[4] = { 71, -8, 3005, 42 };
    int32_t       got[4] = { 0, 0, 0, 0 };
    const h5c_bool_t want_bool[4] = { H5C_TRUE, H5C_TRUE,
                                      H5C_FALSE, H5C_TRUE };
    h5c_bool_t    got_bool[4] = { H5C_FALSE, H5C_FALSE,
                                  H5C_FALSE, H5C_FALSE };
    size_t        count = 0;

    H5C_CHECK(h5c_attr_length(f, "/rank/two", "reference_array", &count));
    H5C_ASSERT_EQ_SIZE(count, 4, "bare-HDF5 attribute length through h5c");
    H5C_CHECK(h5c_read_attr_array(f, "/rank/two", "reference_array", got,
                                  H5C_I32, 4));
    H5C_ASSERT(memcmp(got, want, sizeof want) == 0,
               "bare-HDF5 attribute values differ through h5c");
    H5C_CHECK(h5c_read_attr_array(f, "/rank/two", "reference_bool_array",
                                  got_bool, H5C_BOOL, 4));
    H5C_ASSERT(memcmp(got_bool, want_bool, sizeof want_bool) == 0,
               "integer attribute did not convert to h5c bool values");
}

/* ------------------------------------------------------------------ */
/* phase C: h5c writes                                                 */
/* ------------------------------------------------------------------ */

static void write_with_h5c(void)
{
    h5c_file_t *f = NULL;
    double      d3[K_3D_N];
    const int32_t i32[6] = { 10, -20, 30, -40, 50, -60 };
    const int64_t i64[6] = { 1, -2, 3000000000LL, -4000000000LL, 5, -6 };
    const int16_t i16[6] = { 1, -2, 3, -4, 5, -6 };
    const int8_t  i8[6]  = { 1, -2, 3, -4, 5, -6 };
    const float   f32[6] = { 1.5f, -2.5f, 3.5f, -4.5f, 5.5f, -6.5f };
    const h5c_bool_t flags[6] = { H5C_TRUE, H5C_FALSE, H5C_FALSE,
                                  H5C_TRUE, H5C_TRUE, H5C_FALSE };
    /* Asymmetric per row and per column, so a transposed store shows up. */
    const double u[4] = {   1,   2,   3,   4 };
    const double v[4] = {  10,  20,  30,  40 };
    const double w[4] = { 100, 200, 300, 400 };
    const double *vec[3] = { u, v, w };
    /* Tensor6: six distinct components, XX XY XZ YY YZ ZZ. */
    const double t_xx[2] = { 1, 2 };
    const double t_xy[2] = { 3, 4 };
    const double t_xz[2] = { 5, 6 };
    const double t_yy[2] = { 7, 8 };
    const double t_yz[2] = { 9, 10 };
    const double t_zz[2] = { 11, 12 };
    const double *ten[6] = { t_xx, t_xy, t_xz, t_yy, t_yz, t_zz };
    const size_t  empty_dims[2] = { 0, 3 };
    const double  attr_f64[5] = { -1.25, 8.5, 3.125, -17.0, 42.75 };

    fill_3d(d3);

    H5C_CHECK(h5c_open(H5C_PATH, H5C_TRUNCATE, &f));
    if (f == NULL) {
        return;
    }

    H5C_CHECK(h5c_write_f64(f, "/rank/two", k_2d, 2, k_2d_dims));
    H5C_CHECK(h5c_write_attr_array(f, "/rank/two", "h5c_array", attr_f64,
                                   H5C_F64, 5));
    H5C_CHECK(h5c_write_f64(f, "/rank/three", d3, 3, k_3d_dims));
    H5C_CHECK(h5c_write_f64_scalar(f, "/scalar/r64", 4.0));

    H5C_CHECK(h5c_write_f32(f, "/types/f32", f32, 2, k_2d_dims));
    H5C_CHECK(h5c_write_i8(f, "/types/i8", i8, 2, k_2d_dims));
    H5C_CHECK(h5c_write_i16(f, "/types/i16", i16, 2, k_2d_dims));
    H5C_CHECK(h5c_write_i32(f, "/types/i32", i32, 2, k_2d_dims));
    H5C_CHECK(h5c_write_i64(f, "/types/i64", i64, 2, k_2d_dims));
    H5C_CHECK(h5c_write_bool(f, "/types/bool", flags, 2, k_2d_dims));
    H5C_CHECK(h5c_write_string(f, "/types/text", "hello h5c",
                               H5C_WRITE_DEFAULT));

    H5C_CHECK(h5c_write_interleaved_f64(f, "/fields/velocity", vec, 3, 4));
    H5C_CHECK(h5c_write_interleaved_f64(f, "/fields/stress", ten, 6, 2));

    /* Interleaving by hand must land the same bytes as packing does. */
    {
        const double manual[12] = {
            1, 10, 100,
            2, 20, 200,
            3, 30, 300,
            4, 40, 400
        };
        const size_t dims[2] = { 4, 3 };

        H5C_CHECK(h5c_write_f64(f, "/fields/velocity_manual", manual, 2, dims));
    }

    /* buf may be NULL when the extent is zero. */
    H5C_CHECK(h5c_write_f64(f, "/empty/two", NULL, 2, empty_dims));

    H5C_ASSERT(h5c_file_status(f) == H5C_OK,
               "nothing should have failed while writing " H5C_PATH);
    H5C_CHECK(h5c_close(f));
}

/* ------------------------------------------------------------------ */
/* phase D: the reference reader inspects what h5c wrote               */
/* ------------------------------------------------------------------ */

static void inspect_h5c_bool(hid_t fid)
{
    const unsigned char want_bytes[6] = { 1, 0, 0, 1, 1, 0 };
    hid_t         did, tid, super;
    H5T_class_t   cls;
    int           nmembers, i;
    static const char *const want_names[2]  = { "FALSE", "TRUE" };
    static const int8_t      want_values[2] = { 0, 1 };

    did = H5Dopen2(fid, "/types/bool", H5P_DEFAULT);
    if (did < 0) {
        H5C_FAILF("H5Dopen2 failed for /types/bool");
        return;
    }
    tid = H5Dget_type(did);
    if (tid < 0) {
        H5C_FAILF("H5Dget_type failed for /types/bool");
        H5Dclose(did);
        return;
    }
    cls = H5Tget_class(tid);
    H5C_ASSERT(cls == H5T_ENUM, "bool datatype class %d, want H5T_ENUM (%d)",
               (int)cls, (int)H5T_ENUM);

    if (cls == H5T_ENUM) {
        super = H5Tget_super(tid);
        if (super < 0) {
            H5C_FAILF("H5Tget_super failed for the bool enum");
        } else {
            H5C_ASSERT_EQ_SIZE(H5Tget_size(super), 1,
                               "size of the bool enum's parent type");
            H5Tclose(super);
        }

        nmembers = H5Tget_nmembers(tid);
        H5C_ASSERT(nmembers == 2, "bool enum has %d members, want 2",
                   nmembers);
        for (i = 0; i < nmembers && i < 2; i++) {
            char  *name = H5Tget_member_name(tid, (unsigned)i);
            int8_t value = 0;

            if (name == NULL) {
                H5C_FAILF("H5Tget_member_name failed for member %d", i);
                continue;
            }
            H5C_ASSERT(strcmp(name, want_names[i]) == 0,
                       "bool enum member %d is '%s', want '%s'", i, name,
                       want_names[i]);
            if (H5Tget_member_value(tid, (unsigned)i, &value) < 0) {
                H5C_FAILF("H5Tget_member_value failed for member %d", i);
            } else {
                H5C_ASSERT(value == want_values[i],
                           "bool enum member '%s' has value %d, want %d",
                           name, (int)value, (int)want_values[i]);
            }
            H5free_memory(name);
        }
    }
    H5Tclose(tid);
    H5Dclose(did);

    ref_check_dims(fid, "/types/bool", 2, k_2d_hdims);
    /* One byte per element is the floor: HDF5 elements are byte aligned. */
    ref_check_bytes(fid, "/types/bool", want_bytes, sizeof want_bytes);
}

static void inspect_h5c_string(hid_t fid)
{
    static const char want[] = "hello h5c";
    const size_t      want_len = sizeof want - 1;  /* no NUL on the wire */
    hid_t             did, tid;
    H5T_class_t       cls;
    char              got[sizeof want];

    did = H5Dopen2(fid, "/types/text", H5P_DEFAULT);
    if (did < 0) {
        H5C_FAILF("H5Dopen2 failed for /types/text");
        return;
    }
    tid = H5Dget_type(did);
    if (tid < 0) {
        H5C_FAILF("H5Dget_type failed for /types/text");
        H5Dclose(did);
        return;
    }
    cls = H5Tget_class(tid);
    H5C_ASSERT(cls == H5T_STRING, "string datatype class %d, want %d",
               (int)cls, (int)H5T_STRING);
    if (cls == H5T_STRING) {
        H5C_ASSERT(H5Tis_variable_str(tid) == 0,
                   "the default string form must be fixed length");
        H5C_ASSERT_EQ_SIZE(H5Tget_size(tid), want_len, "stored string length");
        H5C_ASSERT(H5Tget_strpad(tid) == H5T_STR_SPACEPAD,
                   "string padding %d, want H5T_STR_SPACEPAD (%d)",
                   (int)H5Tget_strpad(tid), (int)H5T_STR_SPACEPAD);
        H5C_ASSERT(H5Tget_cset(tid) == H5T_CSET_ASCII,
                   "string charset %d, want H5T_CSET_ASCII (%d)",
                   (int)H5Tget_cset(tid), (int)H5T_CSET_ASCII);

        if (H5Tget_size(tid) == want_len) {
            memset(got, 0, sizeof got);
            if (H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, got) < 0) {
                H5C_FAILF("raw H5Dread failed for /types/text");
            } else {
                H5C_ASSERT(memcmp(got, want, want_len) == 0,
                           "stored string bytes differ from '%s'", want);
            }
        }
    }
    H5Tclose(tid);
    H5Dclose(did);
}

static void inspect_h5c_file(void)
{
    hid_t fid = H5Fopen(H5C_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (fid < 0) {
        H5C_FAILF("H5Fopen failed for " H5C_PATH);
        return;
    }

    /*
     * The core contract: h5c must have stored {3, 2} with the flat byte
     * sequence 1 2 3 4 5 6. Both halves matter -- a reversed implementation
     * that also reordered the buffer would pass the byte check alone.
     */
    ref_check_dims(fid, "/rank/two", 2, k_2d_hdims);
    ref_check_numeric_type(fid, "/rank/two", H5T_FLOAT, 8, 0);
    ref_check_bytes(fid, "/rank/two", k_2d, sizeof k_2d);

    /* Attribute arrays are rank 1 on disk, with the requested extent. */
    {
        const double want[5] = { -1.25, 8.5, 3.125, -17.0, 42.75 };
        double       got[5] = { 0, 0, 0, 0, 0 };
        hsize_t      dims[1] = { 0 };
        hid_t        oid, aid, sid, tid;
        int          rank;

        oid = H5Oopen(fid, "/rank/two", H5P_DEFAULT);
        aid = (oid >= 0) ? H5Aopen(oid, "h5c_array", H5P_DEFAULT)
                         : H5I_INVALID_HID;
        sid = (aid >= 0) ? H5Aget_space(aid) : H5I_INVALID_HID;
        tid = (aid >= 0) ? H5Aget_type(aid) : H5I_INVALID_HID;
        if (oid < 0 || aid < 0 || sid < 0 || tid < 0) {
            H5C_FAILF("cannot inspect h5c array attribute");
        } else {
            rank = H5Sget_simple_extent_ndims(sid);
            H5C_ASSERT(rank == 1, "h5c array attribute rank %d, want 1", rank);
            if (rank == 1 && H5Sget_simple_extent_dims(sid, dims, NULL) >= 0) {
                H5C_ASSERT_EQ_SIZE(dims[0], 5, "h5c array attribute extent");
            } else if (rank == 1) {
                H5C_FAILF("H5Sget_simple_extent_dims failed for h5c_array");
            }
            H5C_ASSERT(H5Tget_class(tid) == H5T_FLOAT &&
                       H5Tget_size(tid) == 8 &&
                       H5Tget_order(tid) == H5T_ORDER_LE,
                       "h5c array attribute is not IEEE F64LE");
            if (H5Aread(aid, H5T_NATIVE_DOUBLE, got) < 0) {
                H5C_FAILF("H5Aread failed for h5c_array");
            } else {
                H5C_ASSERT(memcmp(got, want, sizeof want) == 0,
                           "h5c array attribute values differ in bare HDF5");
            }
        }
        if (tid >= 0) H5Tclose(tid);
        if (sid >= 0) H5Sclose(sid);
        if (aid >= 0) H5Aclose(aid);
        if (oid >= 0) H5Oclose(oid);
    }

    {
        double d3[K_3D_N];

        fill_3d(d3);
        ref_check_dims(fid, "/rank/three", 3, k_3d_hdims);
        ref_check_bytes(fid, "/rank/three", d3, sizeof d3);
    }

    /* A scalar dataset must stay scalar, not become a 1-element array. */
    {
        hid_t did = H5Dopen2(fid, "/scalar/r64", H5P_DEFAULT);

        if (did < 0) {
            H5C_FAILF("H5Dopen2 failed for /scalar/r64");
        } else {
            hid_t sid = H5Dget_space(did);

            H5C_ASSERT(H5Sget_simple_extent_type(sid) == H5S_SCALAR,
                       "/scalar/r64 must be a SCALAR dataspace");
            H5Sclose(sid);
            H5Dclose(did);
        }
    }

    /* Datatypes on disk: class, size and byte order, not just the values. */
    ref_check_numeric_type(fid, "/types/f32", H5T_FLOAT,   4, 0);
    ref_check_numeric_type(fid, "/types/i8",  H5T_INTEGER, 1, 1);
    ref_check_numeric_type(fid, "/types/i16", H5T_INTEGER, 2, 1);
    ref_check_numeric_type(fid, "/types/i32", H5T_INTEGER, 4, 1);
    ref_check_numeric_type(fid, "/types/i64", H5T_INTEGER, 8, 1);

    {
        const int32_t i32[6] = { 10, -20, 30, -40, 50, -60 };
        const int64_t i64[6] = { 1, -2, 3000000000LL, -4000000000LL, 5, -6 };
        const int16_t i16[6] = { 1, -2, 3, -4, 5, -6 };
        const int8_t  i8[6]  = { 1, -2, 3, -4, 5, -6 };
        const float   f32[6] = { 1.5f, -2.5f, 3.5f, -4.5f, 5.5f, -6.5f };

        ref_check_dims(fid, "/types/i32", 2, k_2d_hdims);
        ref_check_bytes(fid, "/types/f32", f32, sizeof f32);
        ref_check_bytes(fid, "/types/i8",  i8,  sizeof i8);
        ref_check_bytes(fid, "/types/i16", i16, sizeof i16);
        ref_check_bytes(fid, "/types/i32", i32, sizeof i32);
        ref_check_bytes(fid, "/types/i64", i64, sizeof i64);
    }

    inspect_h5c_bool(fid);
    inspect_h5c_string(fid);

    /* Interleaved fields: [n, ncomp] with u0 v0 w0 u1 v1 w1 ... */
    {
        const hsize_t vdims[2] = { 4, 3 };
        const double  want[12] = {
            1, 10, 100,
            2, 20, 200,
            3, 30, 300,
            4, 40, 400
        };

        ref_check_dims(fid, "/fields/velocity", 2, vdims);
        ref_check_bytes(fid, "/fields/velocity", want, sizeof want);
        /* Packing by h5c and interleaving by hand must agree byte for byte. */
        ref_check_dims(fid, "/fields/velocity_manual", 2, vdims);
        ref_check_bytes(fid, "/fields/velocity_manual", want, sizeof want);
    }
    {
        /* XDMF3 Tensor6 order XX XY XZ YY YZ ZZ, all six values distinct. */
        const hsize_t tdims[2] = { 2, 6 };
        const double  want[12] = {
            1, 3, 5, 7,  9, 11,
            2, 4, 6, 8, 10, 12
        };

        ref_check_dims(fid, "/fields/stress", 2, tdims);
        ref_check_bytes(fid, "/fields/stress", want, sizeof want);
    }

    /* A zero extent keeps its declared shape and holds no elements. */
    {
        const hsize_t edims[2] = { 0, 3 };
        hid_t         did = H5Dopen2(fid, "/empty/two", H5P_DEFAULT);

        if (did < 0) {
            H5C_FAILF("an empty dataset must still exist as an object");
        } else {
            hid_t sid = H5Dget_space(did);

            H5C_ASSERT(H5Sget_simple_extent_npoints(sid) == 0,
                       "/empty/two must hold no elements");
            H5Sclose(sid);
            H5Dclose(did);
            ref_check_dims(fid, "/empty/two", 2, edims);
        }
    }

    if (H5Fclose(fid) < 0) {
        H5C_FAILF("H5Fclose failed for " H5C_PATH);
    }
}

/* ------------------------------------------------------------------ */
/* phase E: the real h5fortran artifact, when it happens to be there   */
/* ------------------------------------------------------------------ */

/*
 * h5fortran's own test output, read with h5c. This is the only check here
 * that sees bytes a Fortran compiler actually produced, so it is worth
 * having -- but it is a build artifact of a sibling project, so its absence
 * must not fail the suite. h5fortran's tree is never written to.
 */
static void read_h5fortran_artifact(void)
{
    h5c_file_t         *f = NULL;
    h5c_dataset_info_t  info;
    h5c_status_t        st;
    double              got[6];
    const size_t        d22[2] = { 2, 2 };
    const h5c_bool_t    want_lgc[4] = { H5C_TRUE, H5C_FALSE,
                                        H5C_FALSE, H5C_TRUE };
    char               *text = NULL;
    int                 i;

    st = h5c_open(H5FORTRAN_ARTIFACT, H5C_READ, &f);
    if (st != H5C_OK || f == NULL) {
        printf("test_crosslang: h5fortran artifact absent, skipping"
               " (%s)\n", H5FORTRAN_ARTIFACT);
        if (f != NULL) {
            h5c_close(f);
        }
        return;
    }
    printf("test_crosslang: checking against the h5fortran artifact\n");

    /*
     * h5fortran declares real(real64) :: r64_2d(2,3) holding 1..6. The HDF5
     * Fortran library reverses the extents on the way out, so the file holds
     * {3, 2} and the same flat bytes. This is the measurement docs/FORMAT.md
     * is built on.
     */
    H5C_CHECK(h5c_dataset_info(f, "/rank/two", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 2, "h5fortran /rank/two rank");
    H5C_ASSERT_EQ_SIZE(info.dims[0], 3, "h5fortran /rank/two dims[0]");
    H5C_ASSERT_EQ_SIZE(info.dims[1], 2, "h5fortran /rank/two dims[1]");
    H5C_ASSERT(info.type == H5C_F64, "h5fortran /rank/two type %d",
               (int)info.type);

    memset(got, 0, sizeof got);
    H5C_CHECK(h5c_read_f64(f, "/rank/two", got, 2, k_2d_dims));
    for (i = 0; i < 6; i++) {
        H5C_ASSERT(got[i] == k_2d[i],
                   "h5fortran /rank/two flat[%d]: got %g want %g", i, got[i],
                   k_2d[i]);
    }

    /* Reading Fortran-written data with the extents swapped must fail. */
    {
        const size_t swapped[2] = { 2, 3 };

        H5C_CHECK_FAILS(h5c_read_f64(f, "/rank/two", got, 2, swapped),
                        H5C_ERR_SHAPE_MISMATCH);
    }

    /* h5fortran stores `logical` as int32 1/0. See check_int32_logical(). */
    if (h5c_exists(f, "/types/logical")) {
        check_int32_logical(f, "/types/logical", want_lgc, 4, d22);
    }

    /* Fixed-length SPACEPAD written by Fortran, trimmed on the way in. */
    if (h5c_exists(f, "/types/text")) {
        H5C_CHECK(h5c_read_string(f, "/types/text", &text));
        if (text != NULL) {
            H5C_ASSERT(strcmp(text, "hello h5fortran") == 0,
                       "h5fortran /types/text: got '%s'", text);
            h5c_free_string(text);
        }
    }

    /* A string attribute written by h5fortran. */
    if (h5c_attr_exists(f, "/rank/two", "units")) {
        char *units = NULL;

        H5C_CHECK(h5c_read_attr_str(f, "/rank/two", "units", &units));
        if (units != NULL) {
            H5C_ASSERT(strcmp(units, "m/s") == 0,
                       "h5fortran units attribute: got '%s'", units);
            h5c_free_string(units);
        }
    }

    /* The swapped-extent read above is sticky by design. */
    h5c_file_clear_status(f);
    H5C_CHECK(h5c_close(f));
}

/* ------------------------------------------------------------------ */

int main(void)
{
    h5c_file_t *f = NULL;

    H5C_CHECK(h5c_init());

    /*
     * The file image is little-endian by contract. Comparing it against
     * native memory is only meaningful on a little-endian host, so the
     * raw-byte checks stand down elsewhere and the value checks carry on.
     */
    if (!host_is_le()) {
        printf("test_crosslang: big-endian host, raw-byte checks skipped\n");
    }

    /* A: build the reference file with the bare HDF5 C API. */
    build_reference_file();

    /* B: read it back with h5c. */
    H5C_CHECK(h5c_open(REF_PATH, H5C_READ, &f));
    if (f != NULL) {
        read_reference_shapes(f);
        read_reference_types(f);
        read_reference_empty(f);
        read_reference_attr(f);
        H5C_CHECK(h5c_close(f));
    }

    /* C and D: write with h5c, inspect with the bare HDF5 C API. */
    write_with_h5c();
    inspect_h5c_file();

    /* E: optional, only when h5fortran happens to have been built. */
    read_h5fortran_artifact();

    if (g_ref_lcpl != H5I_INVALID_HID) {
        H5Pclose(g_ref_lcpl);
        g_ref_lcpl = H5I_INVALID_HID;
    }
    h5c_finalize();
    return H5C_TEST_SUMMARY("test_crosslang");
}
