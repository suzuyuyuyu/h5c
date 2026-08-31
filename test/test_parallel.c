/*
 * Parallel round-trip tests. Run under mpiexec with 2 ranks by ctest -L quick.
 *
 * Shapes, local extents and values are deliberately ASYMMETRIC: local extents
 * differ between ranks, the non-split dimension is 3, and every element
 * encodes both its rank and its position, so a transposed or mis-offset
 * implementation cannot pass by accident.
 *
 * Every test is bounded by a watchdog alarm: a collective call that fails to
 * agree across ranks must surface as a FAILURE, never as a hung job.
 */
#define _POSIX_C_SOURCE 200809L

#include "h5c_test.h"

#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#include "h5c/h5c_mpi.h"

H5C_TEST_MAIN_STATE;

#define PATH "test_parallel.h5"

/* Wall-clock budget for the whole program; a stall becomes a failed test. */
#define WATCHDOG_SECONDS 120

static int g_me = 0;
static int g_nprocs = 1;

static void on_watchdog(int sig)
{
    (void)sig;
    /* write() is the only safe output here, but the process is dying anyway. */
    fprintf(stderr, "test_parallel: WATCHDOG fired on rank %d "
                    "(likely a collective call that did not agree)\n", g_me);
    fflush(stderr);
    _exit(2);
}

/* Local extents chosen so that no two ranks agree. */
static size_t local_rows(int base)
{
    return (size_t)(base + g_me);
}

/* Expected partition, computed independently of h5c. */
static void expected_partition(size_t nlocal, long long *part)
{
    long long mine = (long long)nlocal;
    int r;

    MPI_Allgather(&mine, 1, MPI_LONG_LONG, part + 1, 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    part[0] = 0;
    for (r = 0; r < g_nprocs; r++) {
        part[r + 1] += part[r];
    }
}

static void check_partition(h5c_file_t *f, const char *path, size_t nlocal)
{
    h5c_dataset_info_t info;
    long long *want;
    int64_t   *got;
    char       dpath[256];
    int        r;

    want = (long long *)calloc((size_t)g_nprocs + 1, sizeof *want);
    got  = (int64_t *)calloc((size_t)g_nprocs + 1, sizeof *got);
    expected_partition(nlocal, want);

    snprintf(dpath, sizeof dpath, "%s/%s", path, H5C_PARTITION_NAME);

    H5C_CHECK(h5c_dataset_info(f, dpath, &info));
    H5C_ASSERT_EQ_SIZE(info.rank, 1, "partition rank");
    H5C_ASSERT_EQ_SIZE(info.dims[0], (size_t)g_nprocs + 1, "partition length");

    H5C_CHECK(h5c_read_i64_1d(f, dpath, got, (size_t)g_nprocs + 1));
    for (r = 0; r <= g_nprocs; r++) {
        H5C_ASSERT((long long)got[r] == want[r],
                   "partition[%d]: got %lld want %lld",
                   r, (long long)got[r], want[r]);
    }
    free(want);
    free(got);
}

/* ------------------------------------------------------------------ */

static void test_1d(h5c_file_t *f)
{
    const size_t nlocal = local_rows(4);   /* 4, 5, 6, ... per rank */
    double      *src, *got;
    long long   *part;
    h5c_dataset_info_t local, global;
    size_t       i;

    src = (double *)malloc(nlocal * sizeof *src);
    got = (double *)malloc(nlocal * sizeof *got);
    part = (long long *)calloc((size_t)g_nprocs + 1, sizeof *part);
    expected_partition(nlocal, part);

    for (i = 0; i < nlocal; i++) {
        src[i] = 1000.0 * g_me + (double)i + 0.25;
        got[i] = -1.0;
    }

    H5C_CHECK(h5c_pwrite(f, "/dist/one", src, H5C_F64, 1, &nlocal,
                         H5C_WRITE_DEFAULT));
    check_partition(f, "/dist/one", nlocal);

    H5C_CHECK(h5c_pdataset_info(f, "/dist/one", &local, &global));
    H5C_ASSERT_EQ_SIZE(local.rank, 1, "1d local rank");
    H5C_ASSERT_EQ_SIZE(local.dims[0], nlocal, "1d local dims[0]");
    H5C_ASSERT_EQ_SIZE(local.count, nlocal, "1d local count");
    H5C_ASSERT_EQ_SIZE(global.rank, 1, "1d global rank");
    H5C_ASSERT_EQ_SIZE(global.dims[0], (size_t)part[g_nprocs],
                       "1d global dims[0]");
    H5C_ASSERT(global.type == H5C_F64, "1d type: got %d", (int)global.type);

    H5C_CHECK(h5c_pread(f, "/dist/one", got, H5C_F64, 1, &nlocal));
    for (i = 0; i < nlocal; i++) {
        H5C_ASSERT(got[i] == src[i], "1d[%lu]: got %g want %g",
                   (unsigned long)i, got[i], src[i]);
    }

    /* The block must sit at this rank's offset, not at the file start. */
    {
        double *whole = (double *)malloc((size_t)part[g_nprocs] * sizeof *whole);
        size_t  total = (size_t)part[g_nprocs];

        H5C_CHECK(h5c_read_f64_1d(f, "/dist/one/data", whole, total));
        for (i = 0; i < nlocal; i++) {
            H5C_ASSERT(whole[(size_t)part[g_me] + i] == src[i],
                       "1d global[%lu]: got %g want %g",
                       (unsigned long)((size_t)part[g_me] + i),
                       whole[(size_t)part[g_me] + i], src[i]);
        }
        free(whole);
    }

    free(src);
    free(got);
    free(part);
}

static void test_2d(h5c_file_t *f)
{
    const size_t nlocal = local_rows(2);   /* 2, 3, 4, ... per rank */
    size_t       ldims[2];
    double      *src, *got;
    h5c_dataset_info_t local, global;
    long long   *part;
    size_t       i, j;

    ldims[0] = nlocal;
    ldims[1] = 3;

    src = (double *)malloc(nlocal * 3 * sizeof *src);
    got = (double *)malloc(nlocal * 3 * sizeof *got);
    part = (long long *)calloc((size_t)g_nprocs + 1, sizeof *part);
    expected_partition(nlocal, part);

    /* Row-major: element (i, j) lives at i*3 + j. Values encode both. */
    for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
            src[i * 3 + j] = 100.0 * g_me + 10.0 * (double)i + (double)j;
        }
    }
    memset(got, 0, nlocal * 3 * sizeof *got);

    H5C_CHECK(h5c_pwrite(f, "/dist/two", src, H5C_F64, 2, ldims,
                         H5C_WRITE_DEFAULT));
    check_partition(f, "/dist/two", nlocal);

    H5C_CHECK(h5c_pdataset_info(f, "/dist/two", &local, &global));
    H5C_ASSERT_EQ_SIZE(local.rank, 2, "2d local rank");
    H5C_ASSERT_EQ_SIZE(local.dims[0], nlocal, "2d local dims[0]");
    H5C_ASSERT_EQ_SIZE(local.dims[1], 3, "2d local dims[1]");
    H5C_ASSERT_EQ_SIZE(local.count, nlocal * 3, "2d local count");
    H5C_ASSERT_EQ_SIZE(global.dims[0], (size_t)part[g_nprocs],
                       "2d global dims[0]");
    H5C_ASSERT_EQ_SIZE(global.dims[1], 3, "2d global dims[1]");

    H5C_CHECK(h5c_pread(f, "/dist/two", got, H5C_F64, 2, ldims));
    for (i = 0; i < nlocal * 3; i++) {
        H5C_ASSERT(got[i] == src[i], "2d flat[%lu]: got %g want %g",
                   (unsigned long)i, got[i], src[i]);
    }

    /* A block placed at the wrong offset would show up here. */
    {
        size_t  total = (size_t)part[g_nprocs];
        size_t  gdims[2];
        double *whole = (double *)malloc(total * 3 * sizeof *whole);

        gdims[0] = total;
        gdims[1] = 3;
        H5C_CHECK(h5c_read_f64(f, "/dist/two/data", whole, 2, gdims));
        for (i = 0; i < nlocal; i++) {
            for (j = 0; j < 3; j++) {
                size_t k = ((size_t)part[g_me] + i) * 3 + j;
                H5C_ASSERT(whole[k] == src[i * 3 + j],
                           "2d global[%lu]: got %g want %g",
                           (unsigned long)k, whole[k], src[i * 3 + j]);
            }
        }
        free(whole);
    }

    free(src);
    free(got);
    free(part);
}

/* An empty rank still joins every collective call. */
static void test_zero_extent(h5c_file_t *f)
{
    const size_t nlocal = (g_me == 0) ? 5 : 0;
    size_t       ldims[2];
    int32_t     *src = NULL, *got = NULL;
    long long   *part;
    h5c_dataset_info_t local, global;
    size_t       i;

    ldims[0] = nlocal;
    ldims[1] = 3;
    part = (long long *)calloc((size_t)g_nprocs + 1, sizeof *part);
    expected_partition(nlocal, part);

    if (nlocal > 0) {
        src = (int32_t *)malloc(nlocal * 3 * sizeof *src);
        got = (int32_t *)malloc(nlocal * 3 * sizeof *got);
        for (i = 0; i < nlocal * 3; i++) {
            src[i] = (int32_t)(7 + 2 * i);
            got[i] = -1;
        }
    }

    H5C_CHECK(h5c_pwrite(f, "/dist/sparse", src, H5C_I32, 2, ldims,
                         H5C_WRITE_DEFAULT));
    check_partition(f, "/dist/sparse", nlocal);

    H5C_CHECK(h5c_pdataset_info(f, "/dist/sparse", &local, &global));
    H5C_ASSERT_EQ_SIZE(local.dims[0], nlocal, "sparse local dims[0]");
    H5C_ASSERT_EQ_SIZE(global.dims[0], 5, "sparse global dims[0]");

    H5C_CHECK(h5c_pread(f, "/dist/sparse", got, H5C_I32, 2, ldims));
    for (i = 0; i < nlocal * 3; i++) {
        H5C_ASSERT(got[i] == src[i], "sparse[%lu]: got %d want %d",
                   (unsigned long)i, (int)got[i], (int)src[i]);
    }

    free(src);
    free(got);
    free(part);
}

/*
 * A disagreement on a non-split dimension must fail on EVERY rank with the
 * same status, without any rank walking into a collective HDF5 call alone.
 * The watchdog turns a deadlock into a failure.
 */
/*
 * Layout accessors: h5c_poffset() and h5c_ppartition().
 *
 * These exist so callers never open "<path>/__partition__" themselves, so the
 * test compares them against a partition computed independently with
 * MPI_Allgather rather than against h5c's own reader.
 */
static void test_layout_accessors(h5c_file_t *f)
{
    /* Asymmetric on purpose: 11, 12, ... rows per rank, non-split extent 2. */
    const size_t nlocal = local_rows(11);
    size_t       ldims[2];
    double      *src;
    long long   *want;
    int64_t     *got;
    size_t       offset = (size_t)-1, mine = (size_t)-1, count = 0;
    size_t       i;
    int          r;

    ldims[0] = nlocal;
    ldims[1] = 2;

    src  = (double *)malloc(nlocal * 2 * sizeof *src);
    want = (long long *)calloc((size_t)g_nprocs + 1, sizeof *want);
    got  = (int64_t *)calloc((size_t)g_nprocs + 1, sizeof *got);
    expected_partition(nlocal, want);

    for (i = 0; i < nlocal * 2; i++) {
        src[i] = 100.0 * g_me + (double)i;
    }
    H5C_CHECK(h5c_pwrite(f, "/layout/field", src, H5C_F64, 2, ldims,
                         H5C_WRITE_DEFAULT));

    /* This rank's slice must match the independently computed boundaries. */
    H5C_CHECK(h5c_poffset(f, "/layout/field", &offset, &mine));
    H5C_ASSERT_EQ_SIZE(offset, (size_t)want[g_me], "poffset offset");
    H5C_ASSERT_EQ_SIZE(mine, nlocal, "poffset local rows");

    /* Either output may be NULL. */
    offset = (size_t)-1;
    H5C_CHECK(h5c_poffset(f, "/layout/field", &offset, NULL));
    H5C_ASSERT_EQ_SIZE(offset, (size_t)want[g_me], "poffset with NULL nlocal");
    H5C_CHECK(h5c_poffset(f, "/layout/field", NULL, &mine));
    H5C_ASSERT_EQ_SIZE(mine, nlocal, "poffset with NULL offset");

    /* Length query first, then the vector itself. */
    H5C_CHECK(h5c_ppartition(f, "/layout/field", NULL, 0, &count));
    H5C_ASSERT_EQ_SIZE(count, (size_t)g_nprocs + 1, "ppartition count");

    H5C_CHECK(h5c_ppartition(f, "/layout/field", got,
                             (size_t)g_nprocs + 1, NULL));
    for (r = 0; r <= g_nprocs; r++) {
        H5C_ASSERT(got[r] == want[r],
                   "ppartition[%d]: got %lld want %lld",
                   r, (long long)got[r], (long long)want[r]);
    }

    /* Too small a buffer is refused, identically on every rank. */
    {
        h5c_status_t st = h5c_ppartition(f, "/layout/field", got,
                                        (size_t)g_nprocs, &count);
        int mn = (int)st, mx = (int)st;
        H5C_ASSERT(st == H5C_ERR_SHAPE_MISMATCH,
                   "short capacity -> %s", h5c_status_string(st));
        MPI_Allreduce(MPI_IN_PLACE, &mn, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        H5C_ASSERT(mn == mx, "short capacity agreed across ranks");
        h5c_file_clear_status(f);
    }

    /* A missing dataset fails, and fails the same way everywhere. */
    {
        h5c_status_t st = h5c_poffset(f, "/layout/absent", &offset, &mine);
        int mn = (int)st, mx = (int)st;
        H5C_ASSERT(st == H5C_ERR_NOT_FOUND,
                   "absent path -> %s", h5c_status_string(st));
        MPI_Allreduce(MPI_IN_PLACE, &mn, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        H5C_ASSERT(mn == mx, "absent path agreed across ranks");
        h5c_file_clear_status(f);
    }

    H5C_ASSERT(h5c_is_parallel(f) == 1, "parallel file reports parallel");

    free(src);
    free(want);
    free(got);
}

/* The accessors must also be right where a rank owns nothing. */
static void test_layout_accessors_empty(h5c_file_t *f)
{
    /* Only the last rank owns rows, so every other offset lands on the end. */
    const size_t nlocal = (g_me == g_nprocs - 1) ? 9 : 0;
    size_t       ldims[2];
    double      *src = NULL;
    long long   *want;
    size_t       offset = (size_t)-1, mine = (size_t)-1;
    size_t       i;

    ldims[0] = nlocal;
    ldims[1] = 2;
    want = (long long *)calloc((size_t)g_nprocs + 1, sizeof *want);
    expected_partition(nlocal, want);

    if (nlocal > 0) {
        src = (double *)malloc(nlocal * 2 * sizeof *src);
        for (i = 0; i < nlocal * 2; i++) {
            src[i] = 500.0 + (double)i;
        }
    }
    /* buf may be NULL exactly when this rank owns no rows. */
    H5C_CHECK(h5c_pwrite(f, "/layout/tail", src, H5C_F64, 2, ldims,
                         H5C_WRITE_DEFAULT));

    H5C_CHECK(h5c_poffset(f, "/layout/tail", &offset, &mine));
    H5C_ASSERT_EQ_SIZE(offset, (size_t)want[g_me], "empty-rank offset");
    H5C_ASSERT_EQ_SIZE(mine, nlocal, "empty-rank local rows");
    if (nlocal == 0) {
        H5C_ASSERT_EQ_SIZE(offset, (size_t)want[g_me + 1],
                           "an empty rank's offset equals its end");
    }

    free(src);
    free(want);
}

static void test_dimension_mismatch(h5c_file_t *f)
{
    size_t  ldims[2];
    double  buf[8];
    h5c_status_t st;
    int     mine, lo, hi;

    if (g_nprocs < 2) {
        return;  /* nothing to disagree about */
    }
    memset(buf, 0, sizeof buf);
    ldims[0] = 2;
    ldims[1] = (size_t)(3 + g_me);   /* differs across ranks, on purpose */

    st = h5c_pwrite(f, "/dist/bad", buf, H5C_F64, 2, ldims,
                    H5C_WRITE_DEFAULT);
    H5C_ASSERT(st == H5C_ERR_SHAPE_MISMATCH,
               "mismatched non-split dim: got %s, expected shape mismatch",
               h5c_status_string(st));

    mine = (int)st;
    MPI_Allreduce(&mine, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mine, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    H5C_ASSERT(lo == hi, "ranks returned different statuses (%d..%d)", lo, hi);

    H5C_ASSERT(!h5c_exists(f, "/dist/bad"),
               "a rejected write must not leave '/dist/bad' behind");
    h5c_file_clear_status(f);

    /* A rank disagreement must also be caught, again on every rank. */
    if (g_me == 0) {
        st = h5c_pwrite(f, "/dist/bad2", buf, H5C_F64, 2, ldims,
                        H5C_WRITE_DEFAULT);
    } else {
        st = h5c_pwrite(f, "/dist/bad2", buf, H5C_F64, 1, ldims,
                        H5C_WRITE_DEFAULT);
    }
    H5C_ASSERT(st == H5C_ERR_SHAPE_MISMATCH,
               "mismatched rank: got %s, expected shape mismatch",
               h5c_status_string(st));
    mine = (int)st;
    MPI_Allreduce(&mine, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mine, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    H5C_ASSERT(lo == hi, "ranks returned different statuses (%d..%d)", lo, hi);
    h5c_file_clear_status(f);
}

/* One rank passes a bad argument; every rank must return together. */
static void test_agreed_invalid_arg(h5c_file_t *f)
{
    const size_t nlocal = 2;
    double       buf[2] = { 1.0, 2.0 };
    h5c_status_t st;
    int          mine, lo, hi;

    if (g_me == 0) {
        st = h5c_pwrite(f, "/dist/argfail", NULL, H5C_F64, 1, &nlocal,
                        H5C_WRITE_DEFAULT);
    } else {
        st = h5c_pwrite(f, "/dist/argfail", buf, H5C_F64, 1, &nlocal,
                        H5C_WRITE_DEFAULT);
    }
    H5C_ASSERT(st == H5C_ERR_INVALID_ARG,
               "one-rank bad argument: got %s, expected invalid argument",
               h5c_status_string(st));
    mine = (int)st;
    MPI_Allreduce(&mine, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mine, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    H5C_ASSERT(lo == hi, "ranks returned different statuses (%d..%d)", lo, hi);
    h5c_file_clear_status(f);
}

static void test_interleaved(h5c_file_t *f)
{
    const size_t n = local_rows(3);   /* 3, 4, 5, ... per rank */
    double      *u, *v, *w, *gu, *gv, *gw;
    const void  *comps[3];
    void        *outs[3];
    h5c_dataset_info_t local, global;
    long long   *part;
    size_t       i;

    u = (double *)malloc(n * sizeof *u);
    v = (double *)malloc(n * sizeof *v);
    w = (double *)malloc(n * sizeof *w);
    gu = (double *)malloc(n * sizeof *gu);
    gv = (double *)malloc(n * sizeof *gv);
    gw = (double *)malloc(n * sizeof *gw);
    part = (long long *)calloc((size_t)g_nprocs + 1, sizeof *part);
    expected_partition(n, part);

    /* Components differ from each other by construction. */
    for (i = 0; i < n; i++) {
        u[i] = 100.0 * g_me + (double)i;
        v[i] = 100.0 * g_me + (double)i + 0.5;
        w[i] = -(100.0 * g_me + (double)i);
        gu[i] = gv[i] = gw[i] = 12345.0;
    }
    comps[0] = u; comps[1] = v; comps[2] = w;
    outs[0] = gu; outs[1] = gv; outs[2] = gw;

    H5C_CHECK(h5c_pwrite_interleaved(f, "/dist/velocity", comps, 3, n,
                                     H5C_F64, H5C_WRITE_DEFAULT));
    check_partition(f, "/dist/velocity", n);

    H5C_CHECK(h5c_pdataset_info(f, "/dist/velocity", &local, &global));
    H5C_ASSERT_EQ_SIZE(local.rank, 2, "interleaved local rank");
    H5C_ASSERT_EQ_SIZE(local.dims[0], n, "interleaved local dims[0]");
    H5C_ASSERT_EQ_SIZE(local.dims[1], 3, "interleaved local dims[1]");
    H5C_ASSERT_EQ_SIZE(global.dims[0], (size_t)part[g_nprocs],
                       "interleaved global dims[0]");
    H5C_ASSERT_EQ_SIZE(global.dims[1], 3, "interleaved global dims[1]");

    H5C_CHECK(h5c_pread_interleaved(f, "/dist/velocity", outs, 3, n,
                                    H5C_F64));
    for (i = 0; i < n; i++) {
        H5C_ASSERT(gu[i] == u[i], "u[%lu]: got %g want %g",
                   (unsigned long)i, gu[i], u[i]);
        H5C_ASSERT(gv[i] == v[i], "v[%lu]: got %g want %g",
                   (unsigned long)i, gv[i], v[i]);
        H5C_ASSERT(gw[i] == w[i], "w[%lu]: got %g want %g",
                   (unsigned long)i, gw[i], w[i]);
    }

    /* The file really is [total, 3] interleaved, not three blocks. */
    {
        size_t  total = (size_t)part[g_nprocs];
        size_t  gdims[2];
        double *whole = (double *)malloc(total * 3 * sizeof *whole);
        size_t  base = (size_t)part[g_me];

        gdims[0] = total;
        gdims[1] = 3;
        H5C_CHECK(h5c_read_f64(f, "/dist/velocity/data", whole, 2, gdims));
        for (i = 0; i < n; i++) {
            H5C_ASSERT(whole[(base + i) * 3 + 0] == u[i],
                       "interleaved u at row %lu", (unsigned long)(base + i));
            H5C_ASSERT(whole[(base + i) * 3 + 1] == v[i],
                       "interleaved v at row %lu", (unsigned long)(base + i));
            H5C_ASSERT(whole[(base + i) * 3 + 2] == w[i],
                       "interleaved w at row %lu", (unsigned long)(base + i));
        }
        free(whole);
    }

    free(u); free(v); free(w);
    free(gu); free(gv); free(gw);
    free(part);
}

/*
 * Tiled interleaved write.
 *
 * The pack limit is small enough to force SEVERAL tiles and the local row
 * counts are so different that the ranks derive DIFFERENT local tile counts
 * (2 versus 4 with two ranks). Only the MPI_Allreduce(MAX) agreement keeps
 * the number of collective calls equal, so a regression here shows up as the
 * watchdog firing rather than as wrong data.
 *
 * The result is compared BYTE FOR BYTE against an untiled reference written
 * from the same components with the default limit.
 */
#define TILED_BASE 200
#define TILED_STEP 320

static void test_interleaved_tiled(h5c_file_t *f)
{
    const size_t n = (size_t)TILED_BASE + (size_t)TILED_STEP * (size_t)g_me;
    double      *u, *v, *w, *gu, *gv, *gw, *ref, *tiled;
    const void  *comps[3];
    void        *outs[3];
    long long   *part;
    size_t       total, base, gdims[2], i, saved;

    u  = (double *)malloc(n * sizeof *u);
    v  = (double *)malloc(n * sizeof *v);
    w  = (double *)malloc(n * sizeof *w);
    gu = (double *)malloc(n * sizeof *gu);
    gv = (double *)malloc(n * sizeof *gv);
    gw = (double *)malloc(n * sizeof *gw);
    part = (long long *)calloc((size_t)g_nprocs + 1, sizeof *part);
    expected_partition(n, part);
    total = (size_t)part[g_nprocs];
    base  = (size_t)part[g_me];
    ref   = (double *)malloc(total * 3 * sizeof *ref);
    tiled = (double *)malloc(total * 3 * sizeof *tiled);

    if (u == NULL || v == NULL || w == NULL || gu == NULL || gv == NULL ||
        gw == NULL || ref == NULL || tiled == NULL) {
        H5C_FAILF("allocation failed in test_interleaved_tiled");
        goto done;
    }

    /* Asymmetric between components, rows and ranks. */
    for (i = 0; i < n; i++) {
        u[i] = 10000.0 * g_me + (double)i + 0.125;
        v[i] = -(10000.0 * g_me + (double)i) - 0.5;
        w[i] = 7.0 * (double)i - 3.0 * g_me;
        gu[i] = gv[i] = gw[i] = 98765.0;
    }
    comps[0] = u; comps[1] = v; comps[2] = w;
    outs[0] = gu; outs[1] = gv; outs[2] = gw;

    saved = h5c_pack_limit();

    /* Reference: one tile per rank at the default limit. */
    H5C_CHECK(h5c_pwrite_interleaved(f, "/dist/tiled_ref", comps, 3, n,
                                     H5C_F64, H5C_WRITE_DEFAULT));

    /* 4096 bytes hold 170 rows of 3 doubles. */
    h5c_set_pack_limit(4096);
    H5C_ASSERT_EQ_SIZE(h5c_pack_limit(), 4096, "pack limit after set");
    H5C_CHECK(h5c_pwrite_interleaved(f, "/dist/tiled", comps, 3, n,
                                     H5C_F64, H5C_WRITE_DEFAULT));
    check_partition(f, "/dist/tiled", n);

    gdims[0] = total;
    gdims[1] = 3;
    H5C_CHECK(h5c_read_f64(f, "/dist/tiled_ref/data", ref, 2, gdims));
    H5C_CHECK(h5c_read_f64(f, "/dist/tiled/data", tiled, 2, gdims));
    for (i = 0; i < total * 3; i++) {
        H5C_ASSERT(tiled[i] == ref[i],
                   "tiled vs untiled flat[%lu]: got %g want %g",
                   (unsigned long)i, tiled[i], ref[i]);
    }
    /* The reference itself must be the interleaved order at this rank's rows. */
    for (i = 0; i < n; i++) {
        H5C_ASSERT(ref[(base + i) * 3 + 0] == u[i] &&
                   ref[(base + i) * 3 + 1] == v[i] &&
                   ref[(base + i) * 3 + 2] == w[i],
                   "tiled row %lu: got %g %g %g want %g %g %g",
                   (unsigned long)(base + i),
                   ref[(base + i) * 3 + 0], ref[(base + i) * 3 + 1],
                   ref[(base + i) * 3 + 2], u[i], v[i], w[i]);
    }

    /* The tiled read must scatter back exactly what was packed. */
    H5C_CHECK(h5c_pread_interleaved(f, "/dist/tiled", outs, 3, n, H5C_F64));
    for (i = 0; i < n; i++) {
        H5C_ASSERT(gu[i] == u[i] && gv[i] == v[i] && gw[i] == w[i],
                   "tiled scatter[%lu]: got %g %g %g want %g %g %g",
                   (unsigned long)i, gu[i], gv[i], gw[i], u[i], v[i], w[i]);
    }

    h5c_set_pack_limit(saved);

done:
    free(u); free(v); free(w);
    free(gu); free(gv); free(gw);
    free(ref); free(tiled);
    free(part);
}

/*
 * A rank with NO local rows inside a tiled interleaved write. Rank 0 owns
 * three tiles' worth of rows, every other rank owns none, so the empty ranks
 * must enter all three collective calls with an empty selection.
 */
static void test_interleaved_tiled_empty(h5c_file_t *f)
{
    const size_t n = (g_me == 0) ? 400 : 0;
    double      *u = NULL, *v = NULL, *gu = NULL, *gv = NULL, *whole = NULL;
    const void  *comps[2];
    void        *outs[2];
    h5c_dataset_info_t local, global;
    size_t       gdims[2], i, saved;

    comps[0] = comps[1] = NULL;
    outs[0] = outs[1] = NULL;
    if (n > 0) {
        u  = (double *)malloc(n * sizeof *u);
        v  = (double *)malloc(n * sizeof *v);
        gu = (double *)malloc(n * sizeof *gu);
        gv = (double *)malloc(n * sizeof *gv);
        if (u == NULL || v == NULL || gu == NULL || gv == NULL) {
            H5C_FAILF("allocation failed in test_interleaved_tiled_empty");
            goto done;
        }
        for (i = 0; i < n; i++) {
            u[i] = (double)i + 0.75;
            v[i] = -2.5 * (double)i - 1.0;
            gu[i] = gv[i] = 54321.0;
        }
        comps[0] = u; comps[1] = v;
        outs[0] = gu; outs[1] = gv;
    }

    saved = h5c_pack_limit();
    /* 4096 bytes hold 256 rows of 2 doubles: rank 0 needs 2 tiles. */
    h5c_set_pack_limit(4096);
    H5C_CHECK(h5c_pwrite_interleaved(f, "/dist/tiled_empty", comps, 2, n,
                                     H5C_F64, H5C_WRITE_DEFAULT));
    check_partition(f, "/dist/tiled_empty", n);

    H5C_CHECK(h5c_pdataset_info(f, "/dist/tiled_empty", &local, &global));
    H5C_ASSERT_EQ_SIZE(local.dims[0], n, "empty-rank local dims[0]");
    H5C_ASSERT_EQ_SIZE(global.dims[0], 400, "empty-rank global dims[0]");
    H5C_ASSERT_EQ_SIZE(global.dims[1], 2, "empty-rank global dims[1]");

    /* Every rank reads the whole dataset and checks the interleaving. */
    whole = (double *)malloc(400 * 2 * sizeof *whole);
    if (whole == NULL) {
        H5C_FAILF("allocation failed in test_interleaved_tiled_empty");
        goto done;
    }
    gdims[0] = 400;
    gdims[1] = 2;
    H5C_CHECK(h5c_read_f64(f, "/dist/tiled_empty/data", whole, 2, gdims));
    for (i = 0; i < 400; i++) {
        double want_u = (double)i + 0.75;
        double want_v = -2.5 * (double)i - 1.0;
        H5C_ASSERT(whole[i * 2 + 0] == want_u && whole[i * 2 + 1] == want_v,
                   "empty-rank row %lu: got %g %g want %g %g",
                   (unsigned long)i, whole[i * 2 + 0], whole[i * 2 + 1],
                   want_u, want_v);
    }

    /* The read side must tolerate the empty local block just the same. */
    H5C_CHECK(h5c_pread_interleaved(f, "/dist/tiled_empty", outs, 2, n,
                                    H5C_F64));
    for (i = 0; i < n; i++) {
        H5C_ASSERT(gu[i] == u[i] && gv[i] == v[i],
                   "empty-rank scatter[%lu]: got %g %g want %g %g",
                   (unsigned long)i, gu[i], gv[i], u[i], v[i]);
    }

    h5c_set_pack_limit(saved);

done:
    free(u); free(v); free(gu); free(gv); free(whole);
}

/*
 * A cross-rank disagreement reached through the interleaved path (different
 * ncomp, hence different dims[1]) must fail on EVERY rank with the same
 * status, even while tiling is active, and leave nothing behind.
 */
static void test_interleaved_mismatch(h5c_file_t *f)
{
    double       a[300], b[300], c[300];
    const void  *comps[3];
    const size_t ncomp = (g_me == 0) ? 3 : 2;
    h5c_status_t st;
    size_t       i, saved;
    int          mine, lo, hi;

    if (g_nprocs < 2) {
        return;  /* nothing to disagree about */
    }
    for (i = 0; i < 300; i++) {
        a[i] = (double)i;
        b[i] = -(double)i;
        c[i] = 0.5 * (double)i;
    }
    comps[0] = a; comps[1] = b; comps[2] = c;

    saved = h5c_pack_limit();
    h5c_set_pack_limit(4096);   /* several tiles, so agreement runs first */
    st = h5c_pwrite_interleaved(f, "/dist/bad_interleaved", comps, ncomp, 300,
                                H5C_F64, H5C_WRITE_DEFAULT);
    H5C_ASSERT(st == H5C_ERR_SHAPE_MISMATCH,
               "mismatched ncomp: got %s, expected shape mismatch",
               h5c_status_string(st));
    mine = (int)st;
    MPI_Allreduce(&mine, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mine, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    H5C_ASSERT(lo == hi, "ranks returned different statuses (%d..%d)", lo, hi);
    H5C_ASSERT(!h5c_exists(f, "/dist/bad_interleaved"),
               "a rejected interleaved write must leave nothing behind");
    h5c_file_clear_status(f);
    h5c_set_pack_limit(saved);
}

static void test_modes(h5c_file_t *f)
{
    H5C_ASSERT(h5c_pis_collective(f) == 1, "collective is the default");
    H5C_ASSERT(h5c_pcomm(f) != MPI_COMM_NULL, "h5c_pcomm returned NULL comm");

    H5C_CHECK(h5c_pset_collective(f, 0));
    H5C_ASSERT(h5c_pis_collective(f) == 0, "independent was not selected");

    /* Independent transfers must still round-trip. */
    {
        const size_t nlocal = local_rows(1);
        double      *src = (double *)malloc(nlocal * sizeof *src);
        double      *got = (double *)malloc(nlocal * sizeof *got);
        size_t       i;

        for (i = 0; i < nlocal; i++) {
            src[i] = 3.5 * (double)(g_me + 1) + (double)i;
            got[i] = 0.0;
        }
        H5C_CHECK(h5c_pwrite(f, "/dist/indep", src, H5C_F64, 1, &nlocal,
                             H5C_WRITE_DEFAULT));
        H5C_CHECK(h5c_pread(f, "/dist/indep", got, H5C_F64, 1, &nlocal));
        for (i = 0; i < nlocal; i++) {
            H5C_ASSERT(got[i] == src[i], "indep[%lu]: got %g want %g",
                       (unsigned long)i, got[i], src[i]);
        }
        free(src);
        free(got);
    }

    H5C_CHECK(h5c_pset_collective(f, 1));
    H5C_ASSERT(h5c_pis_collective(f) == 1, "collective was not restored");
}

/* Reopening read-only exercises the read path against a closed file. */
static void test_reopen_read(void)
{
    h5c_file_t  *f = NULL;
    const size_t nlocal = local_rows(4);
    double      *got;
    size_t       i;

    H5C_CHECK(h5c_popen_comm(PATH, H5C_READ, MPI_COMM_WORLD, MPI_INFO_NULL,
                             &f));
    if (f == NULL) {
        return;
    }
    got = (double *)malloc(nlocal * sizeof *got);
    H5C_CHECK(h5c_pread(f, "/dist/one", got, H5C_F64, 1, &nlocal));
    for (i = 0; i < nlocal; i++) {
        double want = 1000.0 * g_me + (double)i + 0.25;
        H5C_ASSERT(got[i] == want, "reopen 1d[%lu]: got %g want %g",
                   (unsigned long)i, got[i], want);
    }

    /* Writing to a read-only file must fail identically on every rank. */
    {
        h5c_status_t st = h5c_pwrite(f, "/dist/nope", got, H5C_F64, 1,
                                     &nlocal, H5C_WRITE_DEFAULT);
        int mine = (int)st, lo, hi;
        H5C_ASSERT(st == H5C_ERR_STATE, "read-only write: got %s",
                   h5c_status_string(st));
        MPI_Allreduce(&mine, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&mine, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        H5C_ASSERT(lo == hi, "ranks disagreed on the read-only error");
        h5c_file_clear_status(f);
    }

    free(got);
    H5C_CHECK(h5c_close(f));
}

int main(int argc, char **argv)
{
    h5c_file_t *f = NULL;
    int         total = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_me);
    MPI_Comm_size(MPI_COMM_WORLD, &g_nprocs);

    signal(SIGALRM, on_watchdog);
    alarm(WATCHDOG_SECONDS);

    H5C_CHECK(h5c_init());
    H5C_CHECK(h5c_popen(PATH, H5C_TRUNCATE, &f));

    if (f != NULL) {
        test_modes(f);
        test_1d(f);
        test_2d(f);
        test_zero_extent(f);
        test_layout_accessors(f);
        test_layout_accessors_empty(f);
        test_interleaved(f);
        test_interleaved_tiled(f);
        test_interleaved_tiled_empty(f);
        test_interleaved_mismatch(f);
        test_dimension_mismatch(f);
        test_agreed_invalid_arg(f);
        H5C_CHECK(h5c_close(f));
        test_reopen_read();
    }

    h5c_finalize();
    alarm(0);

    /* Any rank's failure must fail the whole test. */
    MPI_Allreduce(&h5c_test_failures, &total, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    h5c_test_failures = total;
    if (g_me != 0) {
        /* Only rank 0 prints the summary; every rank keeps the exit code. */
        MPI_Finalize();
        return (total == 0) ? 0 : 1;
    }
    MPI_Finalize();
    return H5C_TEST_SUMMARY("test_parallel");
}
