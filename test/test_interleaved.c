/*
 * Interleaved multi-component field tests.
 *
 * All data here is ASYMMETRIC: every component uses a different decade, so a
 * transposed, reversed or mis-strided implementation cannot pass by accident.
 * The load-bearing check is the FLAT LAYOUT one: the [n, ncomp] dataset is
 * read back with plain h5c_read_f64 and compared element by element against
 * the expected interleaved byte order.
 */
#include "h5c_test.h"

#include <stdlib.h>
#include <string.h>

H5C_TEST_MAIN_STATE;

#define PATH "test_interleaved.h5"

#define N3 4
#define NCOMP3 3

static const double k_u[N3] = { 1, 2, 3, 4 };
static const double k_v[N3] = { 10, 20, 30, 40 };
static const double k_w[N3] = { 100, 200, 300, 400 };

/* Expected file contents of a {4, 3} dataset built from u, v, w. */
static const double k_flat3[N3 * NCOMP3] = {
    1,  10, 100,
    2,  20, 200,
    3,  30, 300,
    4,  40, 400
};

static void check_flat3(h5c_file_t *f, const char *path)
{
    const size_t dims[2] = { N3, NCOMP3 };
    h5c_dataset_info_t info;
    double got[N3 * NCOMP3];
    int i;

    H5C_CHECK(h5c_dataset_info(f, path, &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 2, "rank");
    H5C_ASSERT_EQ_SIZE(info.dims[0], N3, "dims[0] (rows)");
    H5C_ASSERT_EQ_SIZE(info.dims[1], NCOMP3, "dims[1] (ncomp)");

    memset(got, 0, sizeof got);
    H5C_CHECK(h5c_read_f64(f, path, got, 2, dims));
    for (i = 0; i < N3 * NCOMP3; i++) {
        H5C_ASSERT(got[i] == k_flat3[i],
                   "%s flat[%d]: got %g want %g",
                   path, i, got[i], k_flat3[i]);
    }
}

static void test_vector_roundtrip(h5c_file_t *f)
{
    const void *comps[NCOMP3];
    double gu[N3], gv[N3], gw[N3];
    void *out[NCOMP3];
    int i;

    comps[0] = k_u;
    comps[1] = k_v;
    comps[2] = k_w;
    H5C_CHECK(h5c_write_interleaved(f, "/fields/velocity", comps, NCOMP3, N3,
                                    H5C_F64, H5C_WRITE_DEFAULT));

    /* The check that actually catches ordering bugs. */
    check_flat3(f, "/fields/velocity");

    memset(gu, 0, sizeof gu);
    memset(gv, 0, sizeof gv);
    memset(gw, 0, sizeof gw);
    out[0] = gu;
    out[1] = gv;
    out[2] = gw;
    H5C_CHECK(h5c_read_interleaved(f, "/fields/velocity", out, NCOMP3, N3,
                                   H5C_F64));
    for (i = 0; i < N3; i++) {
        H5C_ASSERT(gu[i] == k_u[i], "u[%d]: got %g want %g", i, gu[i], k_u[i]);
        H5C_ASSERT(gv[i] == k_v[i], "v[%d]: got %g want %g", i, gv[i], k_v[i]);
        H5C_ASSERT(gw[i] == k_w[i], "w[%d]: got %g want %g", i, gw[i], k_w[i]);
    }
}

/* An already-interleaved buffer must produce the identical file contents. */
static void test_matches_plain_write(h5c_file_t *f)
{
    const size_t dims[2] = { N3, NCOMP3 };
    double got[N3 * NCOMP3];
    int i;

    H5C_CHECK(h5c_write_f64(f, "/fields/velocity_plain", k_flat3, 2, dims));
    memset(got, 0, sizeof got);
    H5C_CHECK(h5c_read_f64(f, "/fields/velocity_plain", got, 2, dims));
    for (i = 0; i < N3 * NCOMP3; i++) {
        H5C_ASSERT(got[i] == k_flat3[i],
                   "plain flat[%d]: got %g want %g", i, got[i], k_flat3[i]);
    }
}

static void test_read_component(h5c_file_t *f)
{
    const double *want[NCOMP3];
    double got[N3];
    size_t c;
    int i;

    want[0] = k_u;
    want[1] = k_v;
    want[2] = k_w;

    for (c = 0; c < NCOMP3; c++) {
        memset(got, 0, sizeof got);
        H5C_CHECK(h5c_read_component(f, "/fields/velocity", got, c, N3,
                                     H5C_F64));
        for (i = 0; i < N3; i++) {
            H5C_ASSERT(got[i] == want[c][i],
                       "component %lu [%d]: got %g want %g",
                       (unsigned long)c, i, got[i], want[c][i]);
        }
    }

    /* comp == ncomp is out of range. */
    H5C_CHECK_FAILS(h5c_read_component(f, "/fields/velocity", got, NCOMP3, N3,
                                       H5C_F64),
                    H5C_ERR_INVALID_ARG);
    h5c_file_clear_status(f);
}

/* ncomp = 6: XDMF Tensor6 order XX, XY, XZ, YY, YZ, ZZ, all values distinct. */
#define N6 3
#define NCOMP6 6

static void test_tensor6(h5c_file_t *f)
{
    static const double xx[N6] = { 1, 2, 3 };
    static const double xy[N6] = { 11, 12, 13 };
    static const double xz[N6] = { 21, 22, 23 };
    static const double yy[N6] = { 31, 32, 33 };
    static const double yz[N6] = { 41, 42, 43 };
    static const double zz[N6] = { 51, 52, 53 };
    static const double expect[N6 * NCOMP6] = {
        1, 11, 21, 31, 41, 51,
        2, 12, 22, 32, 42, 52,
        3, 13, 23, 33, 43, 53
    };
    const size_t dims[2] = { N6, NCOMP6 };
    const void *comps[NCOMP6];
    const double *want[NCOMP6];
    double flat[N6 * NCOMP6];
    double got[N6];
    size_t c;
    int i;

    comps[0] = xx; comps[1] = xy; comps[2] = xz;
    comps[3] = yy; comps[4] = yz; comps[5] = zz;
    want[0] = xx; want[1] = xy; want[2] = xz;
    want[3] = yy; want[4] = yz; want[5] = zz;

    H5C_CHECK(h5c_write_interleaved(f, "/fields/stress", comps, NCOMP6, N6,
                                    H5C_F64, H5C_WRITE_DEFAULT));

    memset(flat, 0, sizeof flat);
    H5C_CHECK(h5c_read_f64(f, "/fields/stress", flat, 2, dims));
    for (i = 0; i < N6 * NCOMP6; i++) {
        H5C_ASSERT(flat[i] == expect[i],
                   "tensor6 flat[%d]: got %g want %g", i, flat[i], expect[i]);
    }

    for (c = 0; c < NCOMP6; c++) {
        memset(got, 0, sizeof got);
        H5C_CHECK(h5c_read_component(f, "/fields/stress", got, c, N6, H5C_F64));
        for (i = 0; i < N6; i++) {
            H5C_ASSERT(got[i] == want[c][i],
                       "tensor6 comp %lu [%d]: got %g want %g",
                       (unsigned long)c, i, got[i], want[c][i]);
        }
    }
}

/* ncomp = 1 is the Scalar case and must still yield a {n, 1} dataset. */
static void test_scalar_component(h5c_file_t *f)
{
    const size_t dims[2] = { N3, 1 };
    const void *comps[1];
    double flat[N3];
    double got[N3];
    void *out[1];
    int i;

    comps[0] = k_v;
    H5C_CHECK(h5c_write_interleaved(f, "/fields/pressure", comps, 1, N3,
                                    H5C_F64, H5C_WRITE_DEFAULT));
    H5C_CHECK(h5c_read_f64(f, "/fields/pressure", flat, 2, dims));
    for (i = 0; i < N3; i++) {
        H5C_ASSERT(flat[i] == k_v[i],
                   "scalar flat[%d]: got %g want %g", i, flat[i], k_v[i]);
    }

    memset(got, 0, sizeof got);
    out[0] = got;
    H5C_CHECK(h5c_read_interleaved(f, "/fields/pressure", out, 1, N3, H5C_F64));
    for (i = 0; i < N3; i++) {
        H5C_ASSERT(got[i] == k_v[i],
                   "scalar comp[%d]: got %g want %g", i, got[i], k_v[i]);
    }
}

/* i32 exercises a different element size through the same packing loop. */
static void test_int32(h5c_file_t *f)
{
    static const int32_t a[N3] = { -1, -2, -3, -4 };
    static const int32_t b[N3] = { 7, 8, 9, 10 };
    static const int32_t expect[N3 * 2] = { -1, 7, -2, 8, -3, 9, -4, 10 };
    const size_t dims[2] = { N3, 2 };
    const void *comps[2];
    int32_t flat[N3 * 2];
    int32_t ga[N3], gb[N3];
    void *out[2];
    int i;

    comps[0] = a;
    comps[1] = b;
    H5C_CHECK(h5c_write_interleaved(f, "/fields/pair_i32", comps, 2, N3,
                                    H5C_I32, H5C_WRITE_DEFAULT));
    H5C_CHECK(h5c_read_i32(f, "/fields/pair_i32", flat, 2, dims));
    for (i = 0; i < N3 * 2; i++) {
        H5C_ASSERT(flat[i] == expect[i],
                   "i32 flat[%d]: got %d want %d", i, flat[i], expect[i]);
    }

    out[0] = ga;
    out[1] = gb;
    H5C_CHECK(h5c_read_interleaved(f, "/fields/pair_i32", out, 2, N3, H5C_I32));
    for (i = 0; i < N3; i++) {
        H5C_ASSERT(ga[i] == a[i], "i32 a[%d]: got %d", i, ga[i]);
        H5C_ASSERT(gb[i] == b[i], "i32 b[%d]: got %d", i, gb[i]);
    }
}

/*
 * Tiling must be invisible: with a pack limit small enough to force several
 * tiles, the file contents and the round trip must be identical.
 */
#define NTILE 257

static void test_tiling(h5c_file_t *f)
{
    const size_t dims[2] = { NTILE, NCOMP3 };
    double *p = (double *)malloc(NTILE * sizeof(double));
    double *q = (double *)malloc(NTILE * sizeof(double));
    double *r = (double *)malloc(NTILE * sizeof(double));
    double *ref = (double *)malloc(NTILE * NCOMP3 * sizeof(double));
    double *flat = (double *)malloc(NTILE * NCOMP3 * sizeof(double));
    double *gp = (double *)malloc(NTILE * sizeof(double));
    double *gq = (double *)malloc(NTILE * sizeof(double));
    double *gr = (double *)malloc(NTILE * sizeof(double));
    const void *comps[NCOMP3];
    void *out[NCOMP3];
    size_t saved;
    int i;

    if (p == NULL || q == NULL || r == NULL || ref == NULL || flat == NULL ||
        gp == NULL || gq == NULL || gr == NULL) {
        H5C_FAILF("allocation failed in test_tiling");
        goto done;
    }

    for (i = 0; i < NTILE; i++) {
        p[i] = i + 1;
        q[i] = 1000.0 + i;
        r[i] = -(double)(i + 1) * 3.5;
    }
    comps[0] = p;
    comps[1] = q;
    comps[2] = r;

    /* Reference: one shot, default (large) limit. */
    saved = h5c_pack_limit();
    H5C_ASSERT_EQ_SIZE(saved, (size_t)256u * 1024u * 1024u, "default pack limit");
    H5C_CHECK(h5c_write_interleaved(f, "/tile/ref", comps, NCOMP3, NTILE,
                                    H5C_F64, H5C_WRITE_DEFAULT));
    H5C_CHECK(h5c_read_f64(f, "/tile/ref", ref, 2, dims));

    /* Now force several tiles: 4096 bytes hold 170 rows of 3 doubles. */
    h5c_set_pack_limit(4096);
    H5C_ASSERT_EQ_SIZE(h5c_pack_limit(), 4096, "pack limit after set");
    H5C_CHECK(h5c_write_interleaved(f, "/tile/tiled", comps, NCOMP3, NTILE,
                                    H5C_F64, H5C_WRITE_DEFAULT));
    H5C_CHECK(h5c_read_f64(f, "/tile/tiled", flat, 2, dims));
    for (i = 0; i < NTILE * NCOMP3; i++) {
        H5C_ASSERT(flat[i] == ref[i],
                   "tiled flat[%d]: got %g want %g", i, flat[i], ref[i]);
    }

    /* And the tiled read scatters the same way. */
    out[0] = gp;
    out[1] = gq;
    out[2] = gr;
    H5C_CHECK(h5c_read_interleaved(f, "/tile/tiled", out, NCOMP3, NTILE,
                                   H5C_F64));
    for (i = 0; i < NTILE; i++) {
        H5C_ASSERT(gp[i] == p[i], "tiled p[%d]: got %g want %g", i, gp[i], p[i]);
        H5C_ASSERT(gq[i] == q[i], "tiled q[%d]: got %g want %g", i, gq[i], q[i]);
        H5C_ASSERT(gr[i] == r[i], "tiled r[%d]: got %g want %g", i, gr[i], r[i]);
    }

    /* 0 / 1 are clamped to the minimum rather than rejected. */
    h5c_set_pack_limit(0);
    H5C_ASSERT(h5c_pack_limit() >= 4096, "pack limit clamped to minimum");
    H5C_CHECK(h5c_write_interleaved(f, "/tile/min", comps, NCOMP3, NTILE,
                                    H5C_F64, H5C_WRITE_REPLACE));
    memset(flat, 0, NTILE * NCOMP3 * sizeof(double));
    H5C_CHECK(h5c_read_f64(f, "/tile/min", flat, 2, dims));
    for (i = 0; i < NTILE * NCOMP3; i++) {
        H5C_ASSERT(flat[i] == ref[i],
                   "min-tile flat[%d]: got %g want %g", i, flat[i], ref[i]);
    }

    h5c_set_pack_limit(saved);
    H5C_ASSERT_EQ_SIZE(h5c_pack_limit(), saved, "pack limit restored");

done:
    free(p); free(q); free(r); free(ref); free(flat);
    free(gp); free(gq); free(gr);
}

/*
 * Tiling with a different element size and a tile that does not divide the row
 * count, plus the degenerate zero-row case: both go through the same shared
 * helpers (h5c__tile_rows / h5c__pack_tile / h5c__unpack_tile) that the
 * parallel path uses, so they are worth exercising here too.
 */
#define NTILE_I32 1301

static void test_tiling_edges(h5c_file_t *f)
{
    const size_t dims[2] = { NTILE_I32, NCOMP3 };
    int32_t *a = (int32_t *)malloc(NTILE_I32 * sizeof(int32_t));
    int32_t *b = (int32_t *)malloc(NTILE_I32 * sizeof(int32_t));
    int32_t *c = (int32_t *)malloc(NTILE_I32 * sizeof(int32_t));
    int32_t *flat = (int32_t *)malloc(NTILE_I32 * NCOMP3 * sizeof(int32_t));
    int32_t *ga = (int32_t *)malloc(NTILE_I32 * sizeof(int32_t));
    int32_t *gb = (int32_t *)malloc(NTILE_I32 * sizeof(int32_t));
    int32_t *gc = (int32_t *)malloc(NTILE_I32 * sizeof(int32_t));
    const void *comps[NCOMP3];
    void *out[NCOMP3];
    size_t saved;
    int i;

    if (a == NULL || b == NULL || c == NULL || flat == NULL ||
        ga == NULL || gb == NULL || gc == NULL) {
        H5C_FAILF("allocation failed in test_tiling_edges");
        goto done;
    }

    /* Distinct per component and per row; nothing symmetric. */
    for (i = 0; i < NTILE_I32; i++) {
        a[i] = i + 1;
        b[i] = -(i + 1);
        c[i] = 100000 + 7 * i;
    }
    comps[0] = a;
    comps[1] = b;
    comps[2] = c;

    saved = h5c_pack_limit();
    /* 4096 bytes hold 341 rows of 3 int32: 1301 = 3 * 341 + 278 tiles. */
    h5c_set_pack_limit(4096);
    H5C_CHECK(h5c_write_interleaved(f, "/tile/i32", comps, NCOMP3, NTILE_I32,
                                    H5C_I32, H5C_WRITE_DEFAULT));
    memset(flat, 0, NTILE_I32 * NCOMP3 * sizeof(int32_t));
    H5C_CHECK(h5c_read_i32(f, "/tile/i32", flat, 2, dims));
    for (i = 0; i < NTILE_I32; i++) {
        H5C_ASSERT(flat[i * NCOMP3 + 0] == a[i] &&
                   flat[i * NCOMP3 + 1] == b[i] &&
                   flat[i * NCOMP3 + 2] == c[i],
                   "i32 tiled row %d: got %d %d %d want %d %d %d", i,
                   flat[i * NCOMP3 + 0], flat[i * NCOMP3 + 1],
                   flat[i * NCOMP3 + 2], a[i], b[i], c[i]);
    }

    out[0] = ga;
    out[1] = gb;
    out[2] = gc;
    memset(ga, 0, NTILE_I32 * sizeof(int32_t));
    H5C_CHECK(h5c_read_interleaved(f, "/tile/i32", out, NCOMP3, NTILE_I32,
                                   H5C_I32));
    for (i = 0; i < NTILE_I32; i++) {
        H5C_ASSERT(ga[i] == a[i] && gb[i] == b[i] && gc[i] == c[i],
                   "i32 tiled scatter row %d: got %d %d %d", i,
                   ga[i], gb[i], gc[i]);
    }

    /* Zero rows under a tiny limit must still produce a {0, ncomp} dataset. */
    {
        h5c_dataset_info_t info;

        H5C_CHECK(h5c_write_interleaved(f, "/tile/i32_empty", comps, NCOMP3, 0,
                                        H5C_I32, H5C_WRITE_DEFAULT));
        H5C_CHECK(h5c_dataset_info(f, "/tile/i32_empty", &info));
        H5C_ASSERT_EQ_SIZE(info.dims[0], 0, "rows of empty tiled write");
        H5C_ASSERT_EQ_SIZE(info.dims[1], NCOMP3,
                           "components of empty tiled write");
        H5C_CHECK(h5c_read_interleaved(f, "/tile/i32_empty", out, NCOMP3, 0,
                                       H5C_I32));
    }

    h5c_set_pack_limit(saved);

done:
    free(a); free(b); free(c); free(flat);
    free(ga); free(gb); free(gc);
}

static void test_replace(h5c_file_t *f)
{
    const void *comps[NCOMP3];
    const void *two[2];

    comps[0] = k_u;
    comps[1] = k_v;
    comps[2] = k_w;

    /* Rewriting in place with the same shape is allowed. */
    H5C_CHECK(h5c_write_interleaved(f, "/fields/velocity", comps, NCOMP3, N3,
                                    H5C_F64, H5C_WRITE_DEFAULT));
    check_flat3(f, "/fields/velocity");

    /* A different ncomp changes the shape and needs H5C_WRITE_REPLACE. */
    two[0] = k_u;
    two[1] = k_v;
    H5C_CHECK_FAILS(h5c_write_interleaved(f, "/fields/velocity", two, 2, N3,
                                          H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_SHAPE_MISMATCH);
    h5c_file_clear_status(f);

    H5C_CHECK(h5c_write_interleaved(f, "/fields/velocity", two, 2, N3,
                                    H5C_F64, H5C_WRITE_REPLACE));
    {
        const size_t dims[2] = { N3, 2 };
        const double expect[N3 * 2] = { 1, 10, 2, 20, 3, 30, 4, 40 };
        double got[N3 * 2];
        int i;

        H5C_CHECK(h5c_read_f64(f, "/fields/velocity", got, 2, dims));
        for (i = 0; i < N3 * 2; i++) {
            H5C_ASSERT(got[i] == expect[i],
                       "replaced flat[%d]: got %g want %g", i, got[i], expect[i]);
        }
    }

    /* Put the 3-component field back for later tests. */
    H5C_CHECK(h5c_write_interleaved(f, "/fields/velocity", comps, NCOMP3, N3,
                                    H5C_F64, H5C_WRITE_REPLACE));
    check_flat3(f, "/fields/velocity");
}

static void test_invalid_args(h5c_file_t *f)
{
    const void *comps[NCOMP3];
    const void *bad[NCOMP3];
    void *out[NCOMP3];
    double got[N3];

    comps[0] = k_u;
    comps[1] = k_v;
    comps[2] = k_w;
    bad[0] = k_u;
    bad[1] = NULL;
    bad[2] = k_w;

    H5C_CHECK_FAILS(h5c_write_interleaved(NULL, "/x", comps, NCOMP3, N3,
                                          H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_write_interleaved(f, NULL, comps, NCOMP3, N3,
                                          H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_write_interleaved(f, "/x", NULL, NCOMP3, N3,
                                          H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_write_interleaved(f, "/x", bad, NCOMP3, N3,
                                          H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);
    /* ncomp == 0 is rejected: there is no component to write. */
    H5C_CHECK_FAILS(h5c_write_interleaved(f, "/x", comps, 0, N3,
                                          H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);
    /*
     * n == 0 is LEGAL and yields a {0, ncomp} dataset. Parallel I/O needs an
     * empty local block, so the serial path follows the same rule.
     */
    {
        h5c_dataset_info_t empty_info;
        H5C_CHECK(h5c_write_interleaved(f, "/empty/velocity", comps, NCOMP3, 0,
                                        H5C_F64, H5C_WRITE_DEFAULT));
        H5C_CHECK(h5c_dataset_info(f, "/empty/velocity", &empty_info));
        H5C_ASSERT_EQ_SIZE(empty_info.rank, 2, "rank of empty interleaved");
        H5C_ASSERT_EQ_SIZE(empty_info.dims[0], 0, "rows of empty interleaved");
        H5C_ASSERT_EQ_SIZE(empty_info.dims[1], NCOMP3,
                           "components of empty interleaved");
    }
    /* Strings have no interleaved mapping. */
    H5C_CHECK_FAILS(h5c_write_interleaved(f, "/x", comps, NCOMP3, N3,
                                          H5C_STRING, H5C_WRITE_DEFAULT),
                    H5C_ERR_INVALID_ARG);

    out[0] = got;
    out[1] = NULL;
    out[2] = got;
    H5C_CHECK_FAILS(h5c_read_interleaved(f, "/fields/velocity", out, NCOMP3,
                                         N3, H5C_F64),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_read_interleaved(f, "/fields/velocity", NULL, NCOMP3,
                                         N3, H5C_F64),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_read_component(f, "/fields/velocity", NULL, 0, N3,
                                       H5C_F64),
                    H5C_ERR_INVALID_ARG);
    H5C_CHECK_FAILS(h5c_read_component(f, "/missing", got, 0, N3, H5C_F64),
                    H5C_ERR_NOT_FOUND);
    /* A rank-1 dataset is not an [n, ncomp] field. */
    H5C_CHECK(h5c_write_f64_1d(f, "/plain1d", k_u, N3));
    H5C_CHECK_FAILS(h5c_read_component(f, "/plain1d", got, 0, N3, H5C_F64),
                    H5C_ERR_SHAPE_MISMATCH);
    /* Wrong row count. */
    H5C_CHECK_FAILS(h5c_read_component(f, "/fields/velocity", got, 0, N3 + 1,
                                       H5C_F64),
                    H5C_ERR_SHAPE_MISMATCH);

    h5c_file_clear_status(f);
}

static void test_readonly(void)
{
    h5c_file_t *f = NULL;
    const void *comps[NCOMP3];

    comps[0] = k_u;
    comps[1] = k_v;
    comps[2] = k_w;

    H5C_CHECK(h5c_open(PATH, H5C_READ, &f));
    if (f == NULL) {
        return;
    }
    H5C_CHECK_FAILS(h5c_write_interleaved(f, "/fields/velocity", comps, NCOMP3,
                                          N3, H5C_F64, H5C_WRITE_DEFAULT),
                    H5C_ERR_STATE);
    h5c_file_clear_status(f);

    /* Reading still works on a read-only file. */
    check_flat3(f, "/fields/velocity");
    test_read_component(f);

    h5c_file_clear_status(f);
    H5C_CHECK(h5c_close(f));
}

int main(void)
{
    h5c_file_t *f = NULL;

    H5C_CHECK(h5c_init());

    H5C_CHECK(h5c_open(PATH, H5C_TRUNCATE, &f));
    if (f == NULL) {
        return H5C_TEST_SUMMARY("test_interleaved");
    }

    test_vector_roundtrip(f);
    test_matches_plain_write(f);
    test_read_component(f);
    test_tensor6(f);
    test_scalar_component(f);
    test_int32(f);
    test_tiling(f);
    test_tiling_edges(f);
    test_replace(f);
    test_invalid_args(f);

    H5C_CHECK(h5c_close(f));

    test_readonly();

    h5c_finalize();
    return H5C_TEST_SUMMARY("test_interleaved");
}
