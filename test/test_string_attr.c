/*
 * String datasets and attributes.
 *
 * The fixed-length representation is checked against what h5fortran writes:
 * a scalar H5T_C_S1 of STRSIZE == strlen(value) with H5T_STR_SPACEPAD, so the
 * two libraries can read each other's files.
 */
#include "h5c_test.h"

#include <stdlib.h>
#include <string.h>

H5C_TEST_MAIN_STATE;

#define PATH "test_string_attr.h5"

/* Asserts a round-tripped string and releases the buffer. */
static void expect_string(h5c_file_t *f, const char *path, const char *want)
{
    char *got = NULL;

    H5C_CHECK(h5c_read_string(f, path, &got));
    if (got == NULL) {
        H5C_FAILF("%s: read returned NULL", path);
        return;
    }
    H5C_ASSERT(strcmp(got, want) == 0,
               "%s: got \"%s\", want \"%s\"", path, got, want);
    h5c_free_string(got);
}

static void expect_attr(h5c_file_t *f, const char *obj, const char *name,
                        const char *want)
{
    char *got = NULL;

    H5C_CHECK(h5c_read_attr_str(f, obj, name, &got));
    if (got == NULL) {
        H5C_FAILF("%s/%s: read returned NULL", obj, name);
        return;
    }
    H5C_ASSERT(strcmp(got, want) == 0,
               "%s/%s: got \"%s\", want \"%s\"", obj, name, got, want);
    h5c_free_string(got);
}

/* ------------------------------------------------------------------ */

static void test_fixed_strings(h5c_file_t *f)
{
    h5c_dataset_info_t info;

    H5C_CHECK(h5c_write_string(f, "/text/greeting", "hello h5c",
                               H5C_WRITE_DEFAULT));
    expect_string(f, "/text/greeting", "hello h5c");

    /* Stored as a scalar dataset of string type. */
    H5C_CHECK(h5c_dataset_info(f, "/text/greeting", &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 0, "rank of /text/greeting");
    H5C_ASSERT_EQ_SIZE(info.count, 1, "count of /text/greeting");
    H5C_ASSERT(info.type == H5C_STRING,
               "type of /text/greeting: got %d", (int)info.type);

    /* The empty string is stored as a single space and trims back to "". */
    H5C_CHECK(h5c_write_string(f, "/text/empty", "", H5C_WRITE_DEFAULT));
    expect_string(f, "/text/empty", "");

    /* Trailing spaces cannot survive SPACEPAD; that is the format, not a bug. */
    H5C_CHECK(h5c_write_string(f, "/text/padded", "abc", H5C_WRITE_DEFAULT));
    expect_string(f, "/text/padded", "abc");
}

static void test_grow_and_shrink(h5c_file_t *f)
{
    /* Same length: an in-place rewrite is allowed. */
    H5C_CHECK(h5c_write_string(f, "/text/mutable", "aaaa", H5C_WRITE_DEFAULT));
    H5C_CHECK(h5c_write_string(f, "/text/mutable", "bbbb", H5C_WRITE_DEFAULT));
    expect_string(f, "/text/mutable", "bbbb");

    /* A longer value needs a wider datatype, so it must be requested. */
    H5C_CHECK_FAILS(h5c_write_string(f, "/text/mutable",
                                     "a considerably longer value",
                                     H5C_WRITE_DEFAULT),
                    H5C_ERR_SHAPE_MISMATCH);
    h5c_file_clear_status(f);

    H5C_CHECK(h5c_write_string(f, "/text/mutable",
                               "a considerably longer value",
                               H5C_WRITE_REPLACE));
    expect_string(f, "/text/mutable", "a considerably longer value");

    /* And shrinking again. */
    H5C_CHECK(h5c_write_string(f, "/text/mutable", "x", H5C_WRITE_REPLACE));
    expect_string(f, "/text/mutable", "x");
}

static void test_vlen_strings(h5c_file_t *f)
{
    H5C_CHECK(h5c_write_string_vlen(f, "/text/vlen", "variable length",
                                    H5C_WRITE_DEFAULT));
    expect_string(f, "/text/vlen", "variable length");

    H5C_CHECK(h5c_write_string_vlen(f, "/text/vlen_empty", "",
                                    H5C_WRITE_DEFAULT));
    expect_string(f, "/text/vlen_empty", "");

    /* Without H5C_WRITE_REPLACE an existing vlen string is not overwritten. */
    H5C_CHECK_FAILS(h5c_write_string_vlen(f, "/text/vlen", "other",
                                          H5C_WRITE_DEFAULT),
                    H5C_ERR_EXISTS);
    h5c_file_clear_status(f);

    H5C_CHECK(h5c_write_string_vlen(f, "/text/vlen", "other",
                                    H5C_WRITE_REPLACE));
    expect_string(f, "/text/vlen", "other");
}

static void test_attr_targets(h5c_file_t *f)
{
    const double payload[6] = { 1, 2, 3, 4, 5, 6 };
    const size_t dims[2] = { 3, 2 };

    H5C_CHECK(h5c_write_f64(f, "/rank/two", payload, 2, dims));

    /* A dataset, a group, and the root group. */
    H5C_CHECK(h5c_write_attr_str(f, "/rank/two", "units", "m/s"));
    H5C_CHECK(h5c_write_attr_str(f, "/rank/two", "description",
                                 "asymmetric 3x2 sample"));
    H5C_CHECK(h5c_write_attr_str(f, "/rank", "description", "rank tests"));
    H5C_CHECK(h5c_write_attr_str(f, "/", "creator", "h5c"));

    expect_attr(f, "/rank/two", "units", "m/s");
    expect_attr(f, "/rank/two", "description", "asymmetric 3x2 sample");
    expect_attr(f, "/rank", "description", "rank tests");
    expect_attr(f, "/", "creator", "h5c");

    /* An empty attribute value round-trips too. */
    H5C_CHECK(h5c_write_attr_str(f, "/rank/two", "note", ""));
    expect_attr(f, "/rank/two", "note", "");

    /* Rewriting replaces, including with a longer value. */
    H5C_CHECK(h5c_write_attr_str(f, "/rank/two", "units",
                                 "metres per second"));
    expect_attr(f, "/rank/two", "units", "metres per second");
    H5C_CHECK(h5c_write_attr_str(f, "/rank/two", "units", "m/s"));
    expect_attr(f, "/rank/two", "units", "m/s");
}

static void test_attr_exists(h5c_file_t *f)
{
    H5C_ASSERT(h5c_attr_exists(f, "/rank/two", "units") == 1,
               "attr_exists: /rank/two units should exist");
    H5C_ASSERT(h5c_attr_exists(f, "/rank/two", "missing") == 0,
               "attr_exists: missing attribute reported as present");
    H5C_ASSERT(h5c_attr_exists(f, "/no/such/object", "units") == 0,
               "attr_exists: missing object reported as present");
    H5C_ASSERT(h5c_attr_exists(f, "/", "creator") == 1,
               "attr_exists: root attribute should exist");
    /* Queries never disturb the sticky status. */
    H5C_ASSERT(h5c_file_status(f) == H5C_OK,
               "attr_exists must not record an error");
}

static void test_attr_scalars(h5c_file_t *f)
{
    const double     dt    = 1.25e-3;
    const float      ratio = 0.5f;
    const int32_t    step  = 12345;
    const int64_t    ncell = 4000000000LL;
    const h5c_bool_t done  = H5C_TRUE;

    double     got_dt = 0.0;
    float      got_ratio = 0.0f;
    int32_t    got_step = 0;
    int64_t    got_ncell = 0;
    h5c_bool_t got_done = H5C_FALSE;

    H5C_CHECK(h5c_write_attr_scalar(f, "/rank/two", "dt", &dt, H5C_F64));
    H5C_CHECK(h5c_write_attr_scalar(f, "/rank/two", "ratio", &ratio, H5C_F32));
    H5C_CHECK(h5c_write_attr_scalar(f, "/rank/two", "step", &step, H5C_I32));
    H5C_CHECK(h5c_write_attr_scalar(f, "/rank/two", "ncell", &ncell, H5C_I64));
    H5C_CHECK(h5c_write_attr_scalar(f, "/", "done", &done, H5C_BOOL));

    H5C_CHECK(h5c_read_attr_scalar(f, "/rank/two", "dt", &got_dt, H5C_F64));
    H5C_CHECK(h5c_read_attr_scalar(f, "/rank/two", "ratio", &got_ratio,
                                   H5C_F32));
    H5C_CHECK(h5c_read_attr_scalar(f, "/rank/two", "step", &got_step, H5C_I32));
    H5C_CHECK(h5c_read_attr_scalar(f, "/rank/two", "ncell", &got_ncell,
                                   H5C_I64));
    H5C_CHECK(h5c_read_attr_scalar(f, "/", "done", &got_done, H5C_BOOL));

    H5C_ASSERT(got_dt == dt, "attr dt: got %g want %g", got_dt, dt);
    H5C_ASSERT(got_ratio == ratio, "attr ratio: got %g want %g",
               (double)got_ratio, (double)ratio);
    H5C_ASSERT(got_step == step, "attr step: got %d", (int)got_step);
    H5C_ASSERT(got_ncell == ncell, "attr ncell: got %ld", (long)got_ncell);
    H5C_ASSERT(got_done == H5C_TRUE, "attr done: got %d", (int)got_done);

    /* Overwriting a scalar attribute keeps the newest value. */
    {
        const int32_t step2 = -7;
        int32_t got2 = 0;

        H5C_CHECK(h5c_write_attr_scalar(f, "/rank/two", "step", &step2,
                                        H5C_I32));
        H5C_CHECK(h5c_read_attr_scalar(f, "/rank/two", "step", &got2,
                                       H5C_I32));
        H5C_ASSERT(got2 == step2, "attr step after rewrite: got %d",
                   (int)got2);
    }
}

/*
 * Failure paths. Each of these must leave *out untouched (NULL) so that a
 * caller which frees unconditionally cannot double-free or leak.
 */
static void test_errors(h5c_file_t *f)
{
    char   *s = (char *)0x1; /* poison: a failing read must overwrite it */
    double  v = 0.0;

    H5C_CHECK_FAILS(h5c_read_attr_str(f, "/rank/two", "missing", &s),
                    H5C_ERR_NOT_FOUND);
    H5C_ASSERT(s == NULL, "failed h5c_read_attr_str must set *out to NULL");

    s = (char *)0x1;
    H5C_CHECK_FAILS(h5c_read_attr_str(f, "/no/such/object", "units", &s),
                    H5C_ERR_NOT_FOUND);
    H5C_ASSERT(s == NULL, "failed h5c_read_attr_str must set *out to NULL");

    s = (char *)0x1;
    H5C_CHECK_FAILS(h5c_read_string(f, "/text/missing", &s),
                    H5C_ERR_NOT_FOUND);
    H5C_ASSERT(s == NULL, "failed h5c_read_string must set *out to NULL");

    /* A numeric dataset is not a string. */
    s = (char *)0x1;
    H5C_CHECK_FAILS(h5c_read_string(f, "/rank/two", &s),
                    H5C_ERR_TYPE_MISMATCH);
    H5C_ASSERT(s == NULL, "failed h5c_read_string must set *out to NULL");

    H5C_CHECK_FAILS(h5c_read_attr_scalar(f, "/rank/two", "missing", &v,
                                         H5C_F64),
                    H5C_ERR_NOT_FOUND);

    /* Bad arguments. */
    H5C_CHECK_FAILS(h5c_write_string(f, "/text/null", NULL,
                                     H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_write_attr_str(f, "/rank/two", "", "x"),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_read_string(f, "", &s), H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_write_attr_scalar(f, "/rank/two", "bad", "x",
                                          H5C_STRING),
                    H5C_ERR_INVALID_ARG);

    /* h5c_free_string(NULL) is a no-op, like free(). */
    h5c_free_string(NULL);

    /* Everything above landed on the sticky record by design. */
    h5c_file_clear_status(f);
}

static void test_readonly(void)
{
    h5c_file_t *f = NULL;
    char       *s = NULL;
    double      dt = 0.0;

    H5C_CHECK(h5c_open(PATH, H5C_READ, &f));
    if (f == NULL) {
        return;
    }

    /* Reading still works. */
    H5C_CHECK(h5c_read_string(f, "/text/greeting", &s));
    if (s != NULL) {
        H5C_ASSERT(strcmp(s, "hello h5c") == 0,
                   "read-only reread: got \"%s\"", s);
        h5c_free_string(s);
    }
    H5C_CHECK(h5c_read_attr_scalar(f, "/rank/two", "dt", &dt, H5C_F64));
    H5C_ASSERT(dt == 1.25e-3, "read-only attr dt: got %g", dt);

    /* Writing must be refused before HDF5 is ever touched. */
    H5C_CHECK_FAILS(h5c_write_string(f, "/text/nope", "x",
                                     H5C_WRITE_DEFAULT),
                    H5C_ERR_STATE);
    H5C_CHECK_FAILS(h5c_write_string_vlen(f, "/text/nope", "x",
                                          H5C_WRITE_DEFAULT),
                    H5C_ERR_STATE);
    H5C_CHECK_FAILS(h5c_write_attr_str(f, "/rank/two", "units", "x"),
                    H5C_ERR_STATE);
    H5C_CHECK_FAILS(h5c_write_attr_scalar(f, "/rank/two", "dt", &dt,
                                          H5C_F64),
                    H5C_ERR_STATE);

    h5c_file_clear_status(f);
    H5C_CHECK(h5c_close(f));
}

int main(void)
{
    h5c_file_t *f = NULL;

    H5C_CHECK(h5c_init());

    H5C_CHECK(h5c_open(PATH, H5C_TRUNCATE, &f));
    if (f == NULL) {
        return H5C_TEST_SUMMARY("test_string_attr");
    }

    test_fixed_strings(f);
    test_grow_and_shrink(f);
    test_vlen_strings(f);
    test_attr_targets(f);
    test_attr_exists(f);
    test_attr_scalars(f);
    test_errors(f);

    H5C_CHECK(h5c_close(f));

    test_readonly();

    h5c_finalize();
    return H5C_TEST_SUMMARY("test_string_attr");
}
