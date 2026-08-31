/*
 * Serial round-trip tests for the core h5c API.
 *
 * Shapes and values here are deliberately ASYMMETRIC: a transposed or
 * reversed implementation must not be able to reproduce them by accident.
 */
#include "h5c_test.h"

#include <stdlib.h>
#include <string.h>

H5C_TEST_MAIN_STATE;

#define PATH "test_serial.h5"

/* r64_2d mirrors h5fortran's test: shape (2,3) in Fortran == {3,2} in C. */
static const double k_2d[6] = { 1, 2, 3, 4, 5, 6 };
static const size_t k_2d_dims[2] = { 3, 2 };

static void test_scalar_and_1d(h5c_file_t *f)
{
    const double src1d[3] = { 1.5, -2.5, 3.25 };
    double got1d[3] = { 0, 0, 0 };
    double got_scalar = 0.0;
    int i;

    H5C_CHECK(h5c_write_f64_scalar(f, "/scalar/r64", 42.5));
    H5C_CHECK(h5c_read_f64_scalar(f, "/scalar/r64", &got_scalar));
    H5C_ASSERT(got_scalar == 42.5, "scalar value: got %g", got_scalar);

    H5C_CHECK(h5c_write_f64_1d(f, "/rank/one", src1d, 3));
    H5C_CHECK(h5c_read_f64_1d(f, "/rank/one", got1d, 3));
    for (i = 0; i < 3; i++) {
        H5C_ASSERT(got1d[i] == src1d[i],
                   "1d[%d]: got %g want %g", i, got1d[i], src1d[i]);
    }
}

static void test_multidim_order(h5c_file_t *f)
{
    h5c_dataset_info_t info;
    double got[6];
    int i;

    H5C_CHECK(h5c_write_f64(f, "/rank/two", k_2d, 2, k_2d_dims));

    /* The stored shape must be {3, 2}, not {2, 3}. */
    H5C_CHECK(h5c_dataset_info(f, "/rank/two", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 2, "rank of /rank/two");
    H5C_ASSERT_EQ_SIZE(info.dims[0], 3, "dims[0] of /rank/two");
    H5C_ASSERT_EQ_SIZE(info.dims[1], 2, "dims[1] of /rank/two");
    H5C_ASSERT_EQ_SIZE(info.count, 6, "count of /rank/two");
    H5C_ASSERT(info.type == H5C_F64, "type of /rank/two: got %d", (int)info.type);

    memset(got, 0, sizeof got);
    H5C_CHECK(h5c_read_f64(f, "/rank/two", got, 2, k_2d_dims));
    for (i = 0; i < 6; i++) {
        H5C_ASSERT(got[i] == k_2d[i],
                   "2d flat[%d]: got %g want %g", i, got[i], k_2d[i]);
    }

    /* Reading with the dimensions swapped must be rejected, not reinterpreted. */
    {
        const size_t swapped[2] = { 2, 3 };
        H5C_CHECK_FAILS(h5c_read_f64(f, "/rank/two", got, 2, swapped),
                        H5C_ERR_SHAPE_MISMATCH);
    }

    /* That failure is now on the file's sticky record, by design. */
    h5c_file_clear_status(f);
}

static void test_types(h5c_file_t *f)
{
    const int32_t i32[4] = { 10, -20, 30, -40 };
    const int64_t i64[3] = { 1, -2000000000L, 4000000000L };
    const float   f32[2] = { 0.5f, -1.5f };
    const h5c_bool_t flags[4] = { H5C_TRUE, H5C_FALSE, H5C_FALSE, H5C_TRUE };

    int32_t g32[4];
    int64_t g64[3];
    float   gf32[2];
    h5c_bool_t gflags[4];
    const size_t d22[2] = { 2, 2 };
    int i;

    H5C_CHECK(h5c_write_i32(f, "/types/i32", i32, 2, d22));
    H5C_CHECK(h5c_read_i32(f, "/types/i32", g32, 2, d22));
    for (i = 0; i < 4; i++) {
        H5C_ASSERT(g32[i] == i32[i], "i32[%d]: got %d", i, (int)g32[i]);
    }

    H5C_CHECK(h5c_write_i64_1d(f, "/types/i64", i64, 3));
    H5C_CHECK(h5c_read_i64_1d(f, "/types/i64", g64, 3));
    for (i = 0; i < 3; i++) {
        H5C_ASSERT(g64[i] == i64[i], "i64[%d] mismatch", i);
    }

    H5C_CHECK(h5c_write_f32_1d(f, "/types/f32", f32, 2));
    H5C_CHECK(h5c_read_f32_1d(f, "/types/f32", gf32, 2));
    for (i = 0; i < 2; i++) {
        H5C_ASSERT(gf32[i] == f32[i], "f32[%d] mismatch", i);
    }

    H5C_CHECK(h5c_write_bool(f, "/types/bool", flags, 2, d22));
    H5C_CHECK(h5c_read_bool(f, "/types/bool", gflags, 2, d22));
    for (i = 0; i < 4; i++) {
        H5C_ASSERT(gflags[i] == flags[i], "bool[%d]: got %d", i, (int)gflags[i]);
    }
}

static void test_replace_and_errors(h5c_file_t *f)
{
    const double a[2] = { 1.0, 2.0 };
    const double b[3] = { 7.0, 8.0, 9.0 };
    double got[3];

    H5C_CHECK(h5c_write_f64_1d(f, "/overwrite", a, 2));

    /* Same shape: writing in place is allowed. */
    H5C_CHECK(h5c_write_f64_1d(f, "/overwrite", a, 2));

    /* Different shape without REPLACE is refused. */
    H5C_CHECK_FAILS(h5c_write_f64_1d(f, "/overwrite", b, 3),
                    H5C_ERR_SHAPE_MISMATCH);

    /* With REPLACE the dataset is recreated. */
    {
        const size_t n = 3;
        H5C_CHECK(h5c_write(f, "/overwrite", b, H5C_F64, 1, &n,
                            H5C_WRITE_REPLACE));
        H5C_CHECK(h5c_read_f64_1d(f, "/overwrite", got, 3));
        H5C_ASSERT(got[0] == 7.0 && got[2] == 9.0, "replaced values");
    }

    H5C_CHECK_FAILS(h5c_read_f64_1d(f, "/nope", got, 3), H5C_ERR_NOT_FOUND);
    H5C_ASSERT(h5c_exists(f, "/overwrite") == 1, "exists on a real dataset");
    H5C_ASSERT(h5c_exists(f, "/nope") == 0, "exists on a missing dataset");

    /* Drop the deliberately provoked failures before the caller closes. */
    h5c_file_clear_status(f);
}

static void test_is_parallel(h5c_file_t *f)
{
    /*
     * A serially opened file is never parallel, and the predicate is usable
     * from code that does not include <mpi.h> -- which is the whole point.
     */
    H5C_ASSERT(h5c_is_parallel(f) == 0, "serial file reports not parallel");
    H5C_ASSERT(h5c_is_parallel(NULL) == 0, "NULL handle reports not parallel");
}

static void test_zero_extent(h5c_file_t *f)
{
    /*
     * A zero extent is legal: the dataset exists, carries no elements, and
     * needs no buffer. Parallel I/O depends on this, so serial honours the
     * same rule.
     */
    h5c_dataset_info_t info;
    const size_t empty_1d[1] = { 0 };
    const size_t empty_2d[2] = { 0, 3 };
    double *alloc = NULL;

    H5C_CHECK(h5c_write_f64(f, "/empty/one", NULL, 1, empty_1d));
    H5C_CHECK(h5c_dataset_info(f, "/empty/one", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 1, "rank of /empty/one");
    H5C_ASSERT_EQ_SIZE(info.dims[0], 0, "dims[0] of /empty/one");
    H5C_ASSERT_EQ_SIZE(info.count, 0, "count of /empty/one");
    H5C_CHECK(h5c_read_f64(f, "/empty/one", NULL, 1, empty_1d));

    /* Zero on one axis only; the others keep their extents. */
    H5C_CHECK(h5c_write_f64(f, "/empty/two", NULL, 2, empty_2d));
    H5C_CHECK(h5c_dataset_info(f, "/empty/two", &info));
    H5C_ASSERT_EQ_SIZE(info.dims[0], 0, "dims[0] of /empty/two");
    H5C_ASSERT_EQ_SIZE(info.dims[1], 3, "dims[1] of /empty/two");

    /* read_alloc must still hand back a freeable, non-NULL pointer. */
    H5C_CHECK(h5c_read_alloc(f, "/empty/one", H5C_F64, (void **)&alloc, &info));
    H5C_ASSERT(alloc != NULL, "read_alloc returns non-NULL for an empty dataset");
    h5c_free(alloc);
}

static void test_read_alloc(void)
{
    h5c_file_t *f = NULL;
    h5c_dataset_info_t info;
    double *buf = NULL;
    int i;

    H5C_CHECK(h5c_open(PATH, H5C_READ, &f));
    if (f == NULL) {
        return;
    }
    H5C_CHECK(h5c_read_alloc(f, "/rank/two", H5C_F64, (void **)&buf, &info));
    if (buf != NULL) {
        H5C_ASSERT_EQ_SIZE(info.count, 6, "read_alloc count");
        for (i = 0; i < 6; i++) {
            H5C_ASSERT(buf[i] == k_2d[i], "read_alloc[%d] mismatch", i);
        }
        h5c_free(buf);
    }

    /* A read-only file must refuse writes and remember the failure. */
    H5C_CHECK_FAILS(h5c_write_f64_scalar(f, "/nope", 1.0), H5C_ERR_STATE);
    H5C_ASSERT(h5c_file_status(f) == H5C_ERR_STATE, "sticky status recorded");
    h5c_file_clear_status(f);
    H5C_ASSERT(h5c_file_status(f) == H5C_OK, "sticky status cleared");

    H5C_CHECK(h5c_close(f));
}

static void test_sticky_survives_success(void)
{
    h5c_file_t *f = NULL;
    double v = 0.0;

    H5C_CHECK(h5c_open(PATH, H5C_READ, &f));
    if (f == NULL) {
        return;
    }
    H5C_CHECK_FAILS(h5c_read_f64_scalar(f, "/missing", &v), H5C_ERR_NOT_FOUND);
    /* A later success must NOT clear the earlier failure. */
    H5C_CHECK(h5c_read_f64_scalar(f, "/scalar/r64", &v));
    H5C_ASSERT(h5c_file_status(f) == H5C_ERR_NOT_FOUND,
               "sticky error survived a later success");

    /*
     * h5c_close reports only whether the close succeeded. The sticky error
     * stays visible through h5c_file_status(), which must be read BEFORE
     * closing because the handle is freed here.
     */
    H5C_CHECK(h5c_close(f));
}

int main(void)
{
    h5c_file_t *f = NULL;

    H5C_CHECK(h5c_init());

    H5C_CHECK(h5c_open(PATH, H5C_TRUNCATE, &f));
    if (f == NULL) {
        return H5C_TEST_SUMMARY("test_serial");
    }

    test_scalar_and_1d(f);
    test_multidim_order(f);
    test_types(f);
    test_replace_and_errors(f);
    test_zero_extent(f);
    test_is_parallel(f);

    H5C_CHECK(h5c_close(f));

    test_read_alloc();
    test_sticky_survives_success();

    h5c_finalize();
    return H5C_TEST_SUMMARY("test_serial");
}
