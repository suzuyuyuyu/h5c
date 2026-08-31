/*
 * Dependency-free assertions for the h5c tests.
 *
 * Every test program declares `int h5c_test_failures;` and returns
 * H5C_TEST_SUMMARY(name) from main().
 */
#ifndef H5C_TEST_H
#define H5C_TEST_H

#include <stdio.h>

#include "h5c/h5c.h"

extern int h5c_test_failures;

/* Reports a failure with file/line context and keeps going. */
#define H5C_FAILF(...)                                                        \
    do {                                                                      \
        fprintf(stderr, "  FAIL %s:%d: ", __FILE__, __LINE__);                \
        fprintf(stderr, __VA_ARGS__);                                         \
        fputc('\n', stderr);                                                  \
        h5c_test_failures++;                                                  \
    } while (0)

/* Asserts a boolean condition. */
#define H5C_ASSERT(cond, ...)                                                 \
    do {                                                                      \
        if (!(cond)) {                                                        \
            H5C_FAILF(__VA_ARGS__);                                           \
        }                                                                     \
    } while (0)

/* Asserts that an h5c call returned H5C_OK, printing the recorded message. */
#define H5C_CHECK(expr)                                                       \
    do {                                                                      \
        h5c_status_t h5c_test_st_ = (expr);                                   \
        if (h5c_test_st_ != H5C_OK) {                                         \
            H5C_FAILF("%s -> %s (%s)", #expr,                                 \
                      h5c_status_string(h5c_test_st_),                        \
                      h5c_last_error()->message);                             \
        }                                                                     \
    } while (0)

/* Asserts that an h5c call failed with exactly `want`. */
#define H5C_CHECK_FAILS(expr, want)                                           \
    do {                                                                      \
        h5c_status_t h5c_test_st_ = (expr);                                   \
        if (h5c_test_st_ != (want)) {                                         \
            H5C_FAILF("%s -> %s, expected %s", #expr,                         \
                      h5c_status_string(h5c_test_st_),                        \
                      h5c_status_string(want));                               \
        }                                                                     \
    } while (0)

#define H5C_ASSERT_EQ_SIZE(got, want, what)                                   \
    H5C_ASSERT((size_t)(got) == (size_t)(want), "%s: got %lu, want %lu",       \
               (what), (unsigned long)(got), (unsigned long)(want))

#define H5C_TEST_MAIN_STATE int h5c_test_failures = 0

#define H5C_TEST_SUMMARY(name)                                                \
    (h5c_test_failures == 0                                                   \
         ? (printf("%s: all checks passed\n", (name)), 0)                     \
         : (fprintf(stderr, "%s: %d check(s) failed\n", (name),               \
                    h5c_test_failures), 1))

#endif /* H5C_TEST_H */
