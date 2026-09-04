/*
 * Visualization writer tests. Run under mpiexec (label "mpi"), never on a
 * login node.
 *
 * Everything here is deliberately ASYMMETRIC: per-rank point and cell counts
 * differ, no extent equals another, the three node coordinates live in
 * different decades, and every field component has its own magnitude. A
 * transposed, mis-offset or component-swapped implementation therefore cannot
 * pass by accident.
 *
 * The file is written, closed, reopened and checked against values recomputed
 * from (rank, index) alone. Offsets are recomputed here with MPI_Allgather,
 * which is NOT how h5c_viz derives them (it uses MPI_Exscan), so the two
 * calculations are independent.
 *
 * A watchdog alarm bounds the whole program: a collective call that fails to
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
#include "h5c/h5c_viz.h"

H5C_TEST_MAIN_STATE;

#define PATH "test_pviz.h5"
#define TIME_VALUE 3.75

/* Wall-clock budget for the whole program; a stall becomes a failed test. */
#define WATCHDOG_SECONDS 120

static int g_me = 0;
static int g_nprocs = 1;

static void on_watchdog(int sig)
{
    (void)sig;
    fprintf(stderr, "test_pviz: WATCHDOG fired on rank %d "
                    "(likely a collective call that did not agree)\n", g_me);
    fflush(stderr);
    _exit(2);
}

/* ------------------------------------------------------------------ */
/* the decomposition, recomputed independently of h5c                  */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t np, nc;          /* this rank's counts    */
    size_t poff, coff;      /* this rank's offsets   */
    size_t ptot, ctot;      /* totals over all ranks */
} layout_t;

/* Prefix sums by MPI_Allgather; h5c_viz uses MPI_Exscan, so this is a check. */
static void layout_of(size_t np, size_t nc, layout_t *lay)
{
    long long mine[2];
    long long *all;
    int r;

    all = (long long *)calloc((size_t)g_nprocs * 2, sizeof *all);
    mine[0] = (long long)np;
    mine[1] = (long long)nc;
    MPI_Allgather(mine, 2, MPI_LONG_LONG, all, 2, MPI_LONG_LONG,
                  MPI_COMM_WORLD);

    lay->np = np;
    lay->nc = nc;
    lay->poff = lay->coff = lay->ptot = lay->ctot = 0;
    for (r = 0; r < g_nprocs; r++) {
        if (r < g_me) {
            lay->poff += (size_t)all[2 * r];
            lay->coff += (size_t)all[2 * r + 1];
        }
        lay->ptot += (size_t)all[2 * r];
        lay->ctot += (size_t)all[2 * r + 1];
    }
    free(all);
}

/* ------------------------------------------------------------------ */
/* value generators: pure functions of (rank, index, component)        */
/* ------------------------------------------------------------------ */

/* Three decades apart, so a transposed (3, n) write is unmistakable. */
static double node_c(int comp, int r, size_t i)
{
    switch (comp) {
    case 0:  return 100.0 * r + (double)i + 0.5;
    case 1:  return 100000.0 + 10.0 * (double)i + r;
    default: return -1.0 - (double)r - 0.25 * (double)i;
    }
}

static double pressure_of(int r, size_t i)
{
    return 7.0 * r + 0.125 * (double)i + 1.0;
}

static double velocity_of(size_t comp, int r, size_t i)
{
    return 1000.0 * (double)(comp + 1) + 10.0 * r + (double)i;
}

static double stress_of(size_t comp, int r, size_t i)
{
    return -(double)(comp + 1) - 0.5 * (double)i - 100.0 * r;
}

static int32_t subdomain_of(int r, size_t c)
{
    return (int32_t)(1000 * r + (int)c + 1);
}

static float cellvec_of(size_t comp, int r, size_t c)
{
    return (float)(10 * (int)(comp + 1) + r) + 0.5f * (float)c;
}

static double celltensor_of(size_t comp, int r, size_t c)
{
    return 0.5 + (double)(comp + 1) * 3.0 + 100.0 * r + (double)c;
}

/* Rank-local 0-origin connectivity: cell c, node k. */
static size_t conn_local(size_t c, size_t k, size_t np)
{
    return (c + k) % np;
}

/* ------------------------------------------------------------------ */
/* read-back helpers                                                   */
/* ------------------------------------------------------------------ */

static void check_shape(h5c_file_t *f, const char *path, int rank,
                        size_t d0, size_t d1)
{
    h5c_dataset_info_t info;

    memset(&info, 0, sizeof info);
    H5C_CHECK(h5c_dataset_info(f, path, &info));
    H5C_ASSERT(info.rank == rank, "%s: rank %d, want %d", path, info.rank,
               rank);
    if (info.rank != rank) {
        return;
    }
    H5C_ASSERT_EQ_SIZE(info.dims[0], d0, path);
    if (rank == 2) {
        H5C_ASSERT_EQ_SIZE(info.dims[1], d1, path);
    }
}

static void check_attr_str(h5c_file_t *f, const char *obj, const char *name,
                           const char *want)
{
    char *got = NULL;

    H5C_CHECK(h5c_read_attr_str(f, obj, name, &got));
    if (got == NULL) {
        return;
    }
    H5C_ASSERT(strcmp(got, want) == 0, "%s@%s: got '%s', want '%s'",
               obj, name, got, want);
    h5c_free_string(got);
}

static void check_attr_i32(h5c_file_t *f, const char *obj, const char *name,
                           int32_t want)
{
    int32_t got = -1;

    H5C_CHECK(h5c_read_attr_scalar(f, obj, name, &got, H5C_I32));
    H5C_ASSERT(got == want, "%s@%s: got %d, want %d", obj, name,
               (int)got, (int)want);
}

/*
 * h5c_read() wants the stored rank and extents exactly, so reading an
 * (n, ncomp) dataset needs the 2-D form and a 1-D field the 1-D form. This
 * wraps both shapes and hands back the whole dataset; the caller frees it.
 */
static double *read_f64_rows(h5c_file_t *f, const char *path, size_t rows,
                             size_t ncols)
{
    size_t  dims[2];
    double *buf;
    int     rank = (ncols > 1) ? 2 : 1;

    dims[0] = rows;
    dims[1] = ncols;
    buf = (double *)malloc((rows * ncols > 0 ? rows * ncols : 1) * sizeof *buf);
    if (h5c_read(f, path, buf, H5C_F64, rank, dims) != H5C_OK) {
        H5C_FAILF("cannot read '%s': %s", path, h5c_last_error()->message);
        free(buf);
        return NULL;
    }
    return buf;
}

/* ------------------------------------------------------------------ */
/* mesh 1: unstructured, unequal per-rank point and cell counts         */
/* ------------------------------------------------------------------ */

#define FLUID_NPE 4

static layout_t g_fluid;

static void write_fluid(h5c_viz_t *viz)
{
    h5c_viz_mesh_t mesh;
    double  *nodes, *pressure, *stress, *vel[3], *ctens;
    float   *cvec[3];
    int32_t *conn, *subdom;
    const void *vcomps[3], *ccomps[3];
    size_t   i, c, k, poff = 999, coff = 999;

    layout_of((size_t)(5 + 2 * g_me), (size_t)(2 + g_me), &g_fluid);

    memset(&mesh, 0, sizeof mesh);
    mesh.kind             = H5C_VIZ_UNSTRUCTURED;
    mesh.name             = "fluid";
    mesh.topology         = "Tetrahedron";
    mesh.nodes_per_element = FLUID_NPE;
    mesh.num_points       = g_fluid.np;
    mesh.num_cells        = g_fluid.nc;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    H5C_CHECK(h5c_viz_offsets(viz, &poff, &coff));
    H5C_ASSERT_EQ_SIZE(poff, g_fluid.poff, "fluid point offset");
    H5C_ASSERT_EQ_SIZE(coff, g_fluid.coff, "fluid cell offset");

    /* --- geometry: interleaved nodes, i32 connectivity ------------- */
    nodes = (double *)malloc(g_fluid.np * 3 * sizeof *nodes);
    for (i = 0; i < g_fluid.np; i++) {
        nodes[3 * i + 0] = node_c(0, g_me, i);
        nodes[3 * i + 1] = node_c(1, g_me, i);
        nodes[3 * i + 2] = node_c(2, g_me, i);
    }
    H5C_CHECK(h5c_viz_write_nodes(viz, nodes, H5C_F64));
    free(nodes);

    conn = (int32_t *)malloc(g_fluid.nc * FLUID_NPE * sizeof *conn);
    for (c = 0; c < g_fluid.nc; c++) {
        for (k = 0; k < FLUID_NPE; k++) {
            conn[c * FLUID_NPE + k] =
                (int32_t)conn_local(c, k, g_fluid.np);
        }
    }
    H5C_CHECK(h5c_viz_write_connectivity(viz, conn, H5C_I32));

    /* The caller's buffer must come back untouched: no offset added in place. */
    for (c = 0; c < g_fluid.nc; c++) {
        for (k = 0; k < FLUID_NPE; k++) {
            H5C_ASSERT((size_t)conn[c * FLUID_NPE + k] ==
                           conn_local(c, k, g_fluid.np),
                       "connectivity buffer was modified at [%lu][%lu]",
                       (unsigned long)c, (unsigned long)k);
        }
    }
    free(conn);

    /* --- point data: scalar, vector (comps), 6-tensor (buffer) ----- */
    pressure = (double *)malloc(g_fluid.np * sizeof *pressure);
    for (i = 0; i < g_fluid.np; i++) {
        pressure[i] = pressure_of(g_me, i);
    }
    H5C_CHECK(h5c_viz_write_point_data(viz, "Pressure", pressure, H5C_F64, 1));
    free(pressure);

    for (c = 0; c < 3; c++) {
        vel[c] = (double *)malloc(g_fluid.np * sizeof **vel);
        for (i = 0; i < g_fluid.np; i++) {
            vel[c][i] = velocity_of(c, g_me, i);
        }
        vcomps[c] = vel[c];
    }
    H5C_CHECK(h5c_viz_write_point_data_comps(viz, "Velocity", vcomps,
                                             H5C_F64, 3));
    for (c = 0; c < 3; c++) {
        free(vel[c]);
    }

    stress = (double *)malloc(g_fluid.np * 6 * sizeof *stress);
    for (i = 0; i < g_fluid.np; i++) {
        for (c = 0; c < 6; c++) {
            stress[i * 6 + c] = stress_of(c, g_me, i);
        }
    }
    H5C_CHECK(h5c_viz_write_point_data(viz, "Stress", stress, H5C_F64, 6));
    free(stress);

    /* --- cell data: scalar, vector (comps), 6-tensor (buffer) ------ */
    subdom = (int32_t *)malloc((g_fluid.nc > 0 ? g_fluid.nc : 1)
                               * sizeof *subdom);
    for (c = 0; c < g_fluid.nc; c++) {
        subdom[c] = subdomain_of(g_me, c);
    }
    H5C_CHECK(h5c_viz_write_cell_data(viz, "SubdomainID", subdom, H5C_I32, 1));
    free(subdom);

    for (k = 0; k < 3; k++) {
        cvec[k] = (float *)malloc((g_fluid.nc > 0 ? g_fluid.nc : 1)
                                  * sizeof **cvec);
        for (c = 0; c < g_fluid.nc; c++) {
            cvec[k][c] = cellvec_of(k, g_me, c);
        }
        ccomps[k] = cvec[k];
    }
    H5C_CHECK(h5c_viz_write_cell_data_comps(viz, "CellVelocity", ccomps,
                                            H5C_F32, 3));
    for (k = 0; k < 3; k++) {
        free(cvec[k]);
    }

    ctens = (double *)malloc((g_fluid.nc > 0 ? g_fluid.nc : 1) * 6
                             * sizeof *ctens);
    for (c = 0; c < g_fluid.nc; c++) {
        for (k = 0; k < 6; k++) {
            ctens[c * 6 + k] = celltensor_of(k, g_me, c);
        }
    }
    H5C_CHECK(h5c_viz_write_cell_data(viz, "CellStress", ctens, H5C_F64, 6));
    free(ctens);
}

static void check_fluid(h5c_file_t *f)
{
    double  *nodes, *pressure, *vel, *stress, *ctens;
    float   *cvec;
    int32_t *conn, *subdom;
    size_t   dims[2], i, c, k;

    check_shape(f, "/fluid/geometry/nodes", 2, g_fluid.ptot, 3);
    check_shape(f, "/fluid/geometry/connectivity", 2, g_fluid.ctot, FLUID_NPE);
    check_shape(f, "/fluid/point_data/Pressure", 1, g_fluid.ptot, 0);
    check_shape(f, "/fluid/point_data/Velocity", 2, g_fluid.ptot, 3);
    check_shape(f, "/fluid/point_data/Stress", 2, g_fluid.ptot, 6);
    check_shape(f, "/fluid/cell_data/SubdomainID", 1, g_fluid.ctot, 0);
    check_shape(f, "/fluid/cell_data/CellVelocity", 2, g_fluid.ctot, 3);
    check_shape(f, "/fluid/cell_data/CellStress", 2, g_fluid.ctot, 6);

    check_attr_str(f, "/fluid", "topology_type", "Tetrahedron");
    check_attr_i32(f, "/fluid", "nodes_per_element", FLUID_NPE);
    check_attr_str(f, "/fluid/point_data/Pressure", "attribute_type", "Scalar");
    check_attr_str(f, "/fluid/point_data/Velocity", "attribute_type", "Vector");
    check_attr_str(f, "/fluid/point_data/Stress", "attribute_type", "Tensor6");
    check_attr_str(f, "/fluid/cell_data/SubdomainID", "attribute_type",
                   "Scalar");
    check_attr_str(f, "/fluid/cell_data/CellVelocity", "attribute_type",
                   "Vector");
    check_attr_str(f, "/fluid/cell_data/CellStress", "attribute_type",
                   "Tensor6");

    /* --- nodes: this rank's rows sit at point_offset, x/y/z in order */
    nodes = read_f64_rows(f, "/fluid/geometry/nodes", g_fluid.ptot, 3);
    if (nodes != NULL) {
        for (i = 0; i < g_fluid.np; i++) {
            size_t row = g_fluid.poff + i;
            for (c = 0; c < 3; c++) {
                double want = node_c((int)c, g_me, i);
                H5C_ASSERT(nodes[row * 3 + c] == want,
                           "nodes[%lu][%lu]: got %.17g want %.17g",
                           (unsigned long)row, (unsigned long)c,
                           nodes[row * 3 + c], want);
            }
        }
        free(nodes);
    }

    /*
     * --- connectivity: the writer must have turned rank-local indices into
     * file-global node ids, so every id this rank's cells reference has to
     * fall in [point_offset, point_offset + num_points).
     */
    dims[0] = g_fluid.ctot;
    dims[1] = FLUID_NPE;
    conn = (int32_t *)malloc((g_fluid.ctot * FLUID_NPE > 0
                                  ? g_fluid.ctot * FLUID_NPE : 1)
                             * sizeof *conn);
    H5C_CHECK(h5c_read(f, "/fluid/geometry/connectivity", conn, H5C_I32, 2,
                       dims));
    for (c = 0; c < g_fluid.nc; c++) {
        for (k = 0; k < FLUID_NPE; k++) {
            size_t  row  = g_fluid.coff + c;
            int64_t got  = conn[row * FLUID_NPE + k];
            int64_t want = (int64_t)(g_fluid.poff
                                     + conn_local(c, k, g_fluid.np));

            H5C_ASSERT(got == want,
                       "connectivity[%lu][%lu]: got %lld want %lld",
                       (unsigned long)row, (unsigned long)k,
                       (long long)got, (long long)want);
            H5C_ASSERT(got >= (int64_t)g_fluid.poff &&
                           got < (int64_t)(g_fluid.poff + g_fluid.np),
                       "connectivity[%lu][%lu] = %lld escapes this rank's "
                       "node range [%lu, %lu)",
                       (unsigned long)row, (unsigned long)k, (long long)got,
                       (unsigned long)g_fluid.poff,
                       (unsigned long)(g_fluid.poff + g_fluid.np));
        }
    }
    free(conn);

    /* --- point fields ------------------------------------------------ */
    pressure = read_f64_rows(f, "/fluid/point_data/Pressure", g_fluid.ptot, 1);
    if (pressure != NULL) {
        for (i = 0; i < g_fluid.np; i++) {
            double want = pressure_of(g_me, i);
            H5C_ASSERT(pressure[g_fluid.poff + i] == want,
                       "Pressure[%lu]: got %.17g want %.17g",
                       (unsigned long)(g_fluid.poff + i),
                       pressure[g_fluid.poff + i], want);
        }
        free(pressure);
    }

    vel = read_f64_rows(f, "/fluid/point_data/Velocity", g_fluid.ptot, 3);
    if (vel != NULL) {
        for (i = 0; i < g_fluid.np; i++) {
            for (c = 0; c < 3; c++) {
                double want = velocity_of(c, g_me, i);
                size_t at   = (g_fluid.poff + i) * 3 + c;
                H5C_ASSERT(vel[at] == want,
                           "Velocity[%lu][%lu]: got %.17g want %.17g",
                           (unsigned long)(g_fluid.poff + i),
                           (unsigned long)c, vel[at], want);
            }
        }
        free(vel);
    }

    stress = read_f64_rows(f, "/fluid/point_data/Stress", g_fluid.ptot, 6);
    if (stress != NULL) {
        for (i = 0; i < g_fluid.np; i++) {
            for (c = 0; c < 6; c++) {
                double want = stress_of(c, g_me, i);
                size_t at   = (g_fluid.poff + i) * 6 + c;
                H5C_ASSERT(stress[at] == want,
                           "Stress[%lu][%lu]: got %.17g want %.17g",
                           (unsigned long)(g_fluid.poff + i),
                           (unsigned long)c, stress[at], want);
            }
        }
        free(stress);
    }

    /* --- cell fields ------------------------------------------------- */
    dims[0] = g_fluid.ctot;
    subdom = (int32_t *)malloc((g_fluid.ctot > 0 ? g_fluid.ctot : 1)
                               * sizeof *subdom);
    H5C_CHECK(h5c_read(f, "/fluid/cell_data/SubdomainID", subdom, H5C_I32, 1,
                       dims));
    for (c = 0; c < g_fluid.nc; c++) {
        int32_t want = subdomain_of(g_me, c);
        H5C_ASSERT(subdom[g_fluid.coff + c] == want,
                   "SubdomainID[%lu]: got %d want %d",
                   (unsigned long)(g_fluid.coff + c),
                   (int)subdom[g_fluid.coff + c], (int)want);
    }
    free(subdom);

    dims[0] = g_fluid.ctot;
    dims[1] = 3;
    cvec = (float *)malloc((g_fluid.ctot * 3 > 0 ? g_fluid.ctot * 3 : 1)
                           * sizeof *cvec);
    H5C_CHECK(h5c_read(f, "/fluid/cell_data/CellVelocity", cvec, H5C_F32, 2,
                       dims));
    for (c = 0; c < g_fluid.nc; c++) {
        for (k = 0; k < 3; k++) {
            float  want = cellvec_of(k, g_me, c);
            size_t at   = (g_fluid.coff + c) * 3 + k;
            H5C_ASSERT(cvec[at] == want,
                       "CellVelocity[%lu][%lu]: got %.9g want %.9g",
                       (unsigned long)(g_fluid.coff + c), (unsigned long)k,
                       (double)cvec[at], (double)want);
        }
    }
    free(cvec);

    ctens = read_f64_rows(f, "/fluid/cell_data/CellStress", g_fluid.ctot, 6);
    if (ctens != NULL) {
        for (c = 0; c < g_fluid.nc; c++) {
            for (k = 0; k < 6; k++) {
                double want = celltensor_of(k, g_me, c);
                size_t at   = (g_fluid.coff + c) * 6 + k;
                H5C_ASSERT(ctens[at] == want,
                           "CellStress[%lu][%lu]: got %.17g want %.17g",
                           (unsigned long)(g_fluid.coff + c),
                           (unsigned long)k, ctens[at], want);
            }
        }
        free(ctens);
    }
}

/* ------------------------------------------------------------------ */
/* mesh 2: polydata in the same file as the unstructured one            */
/* ------------------------------------------------------------------ */

static layout_t g_dust;

static void write_dust(h5c_viz_t *viz)
{
    h5c_viz_mesh_t mesh;
    double  *x, *y, *z, *stress;
    const void *xyz[3];
    size_t   i, coff = 999, poff = 999;

    /* Unequal again, and different from the fluid mesh's counts. */
    layout_of((size_t)(3 + 4 * g_me), 0, &g_dust);

    memset(&mesh, 0, sizeof mesh);
    mesh.kind       = H5C_VIZ_POLYDATA;
    mesh.name       = "soil_particles";
    mesh.num_points = g_dust.np;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    H5C_CHECK(h5c_viz_offsets(viz, &poff, &coff));
    H5C_ASSERT_EQ_SIZE(poff, g_dust.poff, "dust point offset");
    H5C_ASSERT_EQ_SIZE(coff, 0, "dust cell offset");

    x = (double *)malloc((g_dust.np > 0 ? g_dust.np : 1) * sizeof *x);
    y = (double *)malloc((g_dust.np > 0 ? g_dust.np : 1) * sizeof *y);
    z = (double *)malloc((g_dust.np > 0 ? g_dust.np : 1) * sizeof *z);
    for (i = 0; i < g_dust.np; i++) {
        x[i] = node_c(0, g_me, i);
        y[i] = node_c(1, g_me, i);
        z[i] = node_c(2, g_me, i);
    }
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    H5C_CHECK(h5c_viz_write_nodes_comps(viz, xyz, H5C_F64));
    free(x);
    free(y);
    free(z);

    stress = (double *)malloc((g_dust.np > 0 ? g_dust.np : 1) * sizeof *stress);
    for (i = 0; i < g_dust.np; i++) {
        stress[i] = pressure_of(g_me, i);
    }
    H5C_CHECK(h5c_viz_write_point_data(viz, "EquivalentStress", stress,
                                       H5C_F64, 1));
    free(stress);
}

static void check_dust(h5c_file_t *f)
{
    double *nodes;
    size_t  i, c;

    check_shape(f, "/soil_particles/geometry/nodes", 2, g_dust.ptot, 3);
    check_shape(f, "/soil_particles/point_data/EquivalentStress", 1,
                g_dust.ptot, 0);
    check_attr_str(f, "/soil_particles", "topology_type", "Polyvertex");
    check_attr_i32(f, "/soil_particles", "nodes_per_element", 1);
    H5C_ASSERT(!h5c_exists(f, "/soil_particles/geometry/connectivity"),
               "a point cloud must have no connectivity dataset");
    H5C_ASSERT(!h5c_exists(f, "/soil_particles/cell_data"),
               "a point cloud must have no cell_data group");

    /* The separate-coordinate form must interleave, not concatenate. */
    nodes = read_f64_rows(f, "/soil_particles/geometry/nodes", g_dust.ptot, 3);
    if (nodes != NULL) {
        for (i = 0; i < g_dust.np; i++) {
            for (c = 0; c < 3; c++) {
                double want = node_c((int)c, g_me, i);
                size_t at   = (g_dust.poff + i) * 3 + c;
                H5C_ASSERT(nodes[at] == want,
                           "dust nodes[%lu][%lu]: got %.17g want %.17g",
                           (unsigned long)(g_dust.poff + i),
                           (unsigned long)c, nodes[at], want);
            }
        }
        free(nodes);
    }
}

/* ------------------------------------------------------------------ */
/* mesh 3: rank 0 owns zero points and zero cells                      */
/* ------------------------------------------------------------------ */

#define SPARSE_NPE 3

static layout_t g_sparse;

static void write_sparse(h5c_viz_t *viz)
{
    h5c_viz_mesh_t mesh;
    double  *nodes;
    int64_t *conn;
    double  *field;
    int32_t *cfield;
    size_t   i, c, k;

    /* Rank 0 contributes nothing at all and must still be collective. */
    layout_of(g_me == 0 ? 0 : (size_t)(4 + g_me),
              g_me == 0 ? 0 : (size_t)(1 + g_me), &g_sparse);

    memset(&mesh, 0, sizeof mesh);
    mesh.kind             = H5C_VIZ_UNSTRUCTURED;
    mesh.name             = "sparse";
    mesh.topology         = "Triangle";
    mesh.nodes_per_element = SPARSE_NPE;
    mesh.num_points       = g_sparse.np;
    mesh.num_cells        = g_sparse.nc;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    /* An empty rank may pass NULL buffers; that is guaranteed, not tolerated. */
    if (g_sparse.np == 0) {
        H5C_CHECK(h5c_viz_write_nodes(viz, NULL, H5C_F64));
    } else {
        nodes = (double *)malloc(g_sparse.np * 3 * sizeof *nodes);
        for (i = 0; i < g_sparse.np; i++) {
            for (c = 0; c < 3; c++) {
                nodes[i * 3 + c] = node_c((int)c, g_me, i);
            }
        }
        H5C_CHECK(h5c_viz_write_nodes(viz, nodes, H5C_F64));
        free(nodes);
    }

    if (g_sparse.nc == 0) {
        H5C_CHECK(h5c_viz_write_connectivity(viz, NULL, H5C_I64));
    } else {
        conn = (int64_t *)malloc(g_sparse.nc * SPARSE_NPE * sizeof *conn);
        for (c = 0; c < g_sparse.nc; c++) {
            for (k = 0; k < SPARSE_NPE; k++) {
                conn[c * SPARSE_NPE + k] =
                    (int64_t)conn_local(c, k, g_sparse.np);
            }
        }
        H5C_CHECK(h5c_viz_write_connectivity(viz, conn, H5C_I64));
        free(conn);
    }

    field = (double *)malloc((g_sparse.np > 0 ? g_sparse.np : 1)
                             * sizeof *field);
    for (i = 0; i < g_sparse.np; i++) {
        field[i] = pressure_of(g_me, i);
    }
    H5C_CHECK(h5c_viz_write_point_data(viz, "Pressure",
                                       g_sparse.np > 0 ? field : NULL,
                                       H5C_F64, 1));
    free(field);

    cfield = (int32_t *)malloc((g_sparse.nc > 0 ? g_sparse.nc : 1)
                               * sizeof *cfield);
    for (c = 0; c < g_sparse.nc; c++) {
        cfield[c] = subdomain_of(g_me, c);
    }
    H5C_CHECK(h5c_viz_write_cell_data(viz, "SubdomainID",
                                      g_sparse.nc > 0 ? cfield : NULL,
                                      H5C_I32, 1));
    free(cfield);
}

static void check_sparse(h5c_file_t *f)
{
    double  *nodes, *field;
    int64_t *conn;
    size_t   dims[2], i, c, k;

    check_shape(f, "/sparse/geometry/nodes", 2, g_sparse.ptot, 3);
    check_shape(f, "/sparse/geometry/connectivity", 2, g_sparse.ctot,
                SPARSE_NPE);
    check_shape(f, "/sparse/point_data/Pressure", 1, g_sparse.ptot, 0);
    check_shape(f, "/sparse/cell_data/SubdomainID", 1, g_sparse.ctot, 0);
    check_attr_str(f, "/sparse", "topology_type", "Triangle");
    check_attr_i32(f, "/sparse", "nodes_per_element", SPARSE_NPE);

    nodes = read_f64_rows(f, "/sparse/geometry/nodes", g_sparse.ptot, 3);
    if (nodes != NULL) {
        for (i = 0; i < g_sparse.np; i++) {
            for (c = 0; c < 3; c++) {
                double want = node_c((int)c, g_me, i);
                size_t at   = (g_sparse.poff + i) * 3 + c;
                H5C_ASSERT(nodes[at] == want,
                           "sparse nodes[%lu][%lu]: got %.17g want %.17g",
                           (unsigned long)(g_sparse.poff + i),
                           (unsigned long)c, nodes[at], want);
            }
        }
        free(nodes);
    }

    dims[0] = g_sparse.ctot;
    dims[1] = SPARSE_NPE;
    conn = (int64_t *)malloc((g_sparse.ctot * SPARSE_NPE > 0
                                  ? g_sparse.ctot * SPARSE_NPE : 1)
                             * sizeof *conn);
    H5C_CHECK(h5c_read(f, "/sparse/geometry/connectivity", conn, H5C_I64, 2,
                       dims));
    for (c = 0; c < g_sparse.nc; c++) {
        for (k = 0; k < SPARSE_NPE; k++) {
            size_t  at   = (g_sparse.coff + c) * SPARSE_NPE + k;
            int64_t want = (int64_t)(g_sparse.poff
                                     + conn_local(c, k, g_sparse.np));
            H5C_ASSERT(conn[at] == want,
                       "sparse connectivity[%lu][%lu]: got %lld want %lld",
                       (unsigned long)(g_sparse.coff + c), (unsigned long)k,
                       (long long)conn[at], (long long)want);
        }
    }
    free(conn);

    field = read_f64_rows(f, "/sparse/point_data/Pressure", g_sparse.ptot, 1);
    if (field != NULL) {
        for (i = 0; i < g_sparse.np; i++) {
            double want = pressure_of(g_me, i);
            H5C_ASSERT(field[g_sparse.poff + i] == want,
                       "sparse Pressure[%lu]: got %.17g want %.17g",
                       (unsigned long)(g_sparse.poff + i),
                       field[g_sparse.poff + i], want);
        }
        free(field);
    }
}

/* ------------------------------------------------------------------ */
/* connectivity in every supported integer width                       */
/* ------------------------------------------------------------------ */

/*
 * i8 has to survive the offset addition, so the meshes here stay small
 * enough that a global node id fits in a signed byte at any rank count that
 * a test job would use.
 */
#define NARROW_NPE 3

typedef struct {
    const char *name;
    h5c_type_t  type;
    layout_t    lay;
} narrow_case_t;

static narrow_case_t g_narrow[4] = {
    { "conn_i8",  H5C_I8,  { 0, 0, 0, 0, 0, 0 } },
    { "conn_i16", H5C_I16, { 0, 0, 0, 0, 0, 0 } },
    { "conn_i32", H5C_I32, { 0, 0, 0, 0, 0, 0 } },
    { "conn_i64", H5C_I64, { 0, 0, 0, 0, 0, 0 } }
};

static void write_narrow(h5c_viz_t *viz, narrow_case_t *nc)
{
    h5c_viz_mesh_t mesh;
    double  *nodes;
    int8_t  *c8;
    int16_t *c16;
    int32_t *c32;
    int64_t *c64;
    void    *conn;
    size_t   n, i, c, k;

    layout_of((size_t)(3 + g_me), (size_t)(1 + 2 * g_me), &nc->lay);

    memset(&mesh, 0, sizeof mesh);
    mesh.kind             = H5C_VIZ_UNSTRUCTURED;
    mesh.name             = nc->name;
    mesh.topology         = "Triangle";
    mesh.nodes_per_element = NARROW_NPE;
    mesh.num_points       = nc->lay.np;
    mesh.num_cells        = nc->lay.nc;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    nodes = (double *)malloc(nc->lay.np * 3 * sizeof *nodes);
    for (i = 0; i < nc->lay.np; i++) {
        for (c = 0; c < 3; c++) {
            nodes[i * 3 + c] = node_c((int)c, g_me, i);
        }
    }
    H5C_CHECK(h5c_viz_write_nodes(viz, nodes, H5C_F64));
    free(nodes);

    n    = nc->lay.nc * NARROW_NPE;
    conn = malloc((n > 0 ? n : 1) * h5c_type_size(nc->type));
    c8   = (int8_t  *)conn;
    c16  = (int16_t *)conn;
    c32  = (int32_t *)conn;
    c64  = (int64_t *)conn;
    for (c = 0; c < nc->lay.nc; c++) {
        for (k = 0; k < NARROW_NPE; k++) {
            size_t at = c * NARROW_NPE + k;
            size_t v  = conn_local(c, k, nc->lay.np);
            switch (nc->type) {
            case H5C_I8:  c8[at]  = (int8_t)v;  break;
            case H5C_I16: c16[at] = (int16_t)v; break;
            case H5C_I32: c32[at] = (int32_t)v; break;
            default:      c64[at] = (int64_t)v; break;
            }
        }
    }
    H5C_CHECK(h5c_viz_write_connectivity(viz, conn, nc->type));
    free(conn);
}

static void check_narrow(h5c_file_t *f, const narrow_case_t *nc)
{
    h5c_dataset_info_t info;
    char    path[128];
    int64_t *conn;
    size_t   dims[2], c, k;

    snprintf(path, sizeof path, "/%s/geometry/connectivity", nc->name);
    check_shape(f, path, 2, nc->lay.ctot, NARROW_NPE);

    memset(&info, 0, sizeof info);
    H5C_CHECK(h5c_dataset_info(f, path, &info));
    H5C_ASSERT(info.type == nc->type,
               "%s: stored type %d, want %d", path, (int)info.type,
               (int)nc->type);

    /* HDF5 widens on read, so one int64 check covers every stored width. */
    dims[0] = nc->lay.ctot;
    dims[1] = NARROW_NPE;
    conn = (int64_t *)malloc((nc->lay.ctot * NARROW_NPE > 0
                                  ? nc->lay.ctot * NARROW_NPE : 1)
                             * sizeof *conn);
    H5C_CHECK(h5c_read(f, path, conn, H5C_I64, 2, dims));
    for (c = 0; c < nc->lay.nc; c++) {
        for (k = 0; k < NARROW_NPE; k++) {
            size_t  at   = (nc->lay.coff + c) * NARROW_NPE + k;
            int64_t want = (int64_t)(nc->lay.poff
                                     + conn_local(c, k, nc->lay.np));
            H5C_ASSERT(conn[at] == want,
                       "%s[%lu][%lu]: got %lld want %lld", path,
                       (unsigned long)(nc->lay.coff + c), (unsigned long)k,
                       (long long)conn[at], (long long)want);
            H5C_ASSERT(conn[at] >= (int64_t)nc->lay.poff &&
                           conn[at] < (int64_t)(nc->lay.poff + nc->lay.np),
                       "%s[%lu][%lu] = %lld escapes this rank's node range",
                       path, (unsigned long)(nc->lay.coff + c),
                       (unsigned long)k, (long long)conn[at]);
        }
    }
    free(conn);
}

/* ------------------------------------------------------------------ */
/* rejected arguments, agreed across ranks                             */
/* ------------------------------------------------------------------ */

/* Asserts that every rank returned the same status. */
static void check_agreed(h5c_status_t st, const char *what)
{
    int mine = (int)st, lo = 0, hi = 0;

    MPI_Allreduce(&mine, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mine, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    H5C_ASSERT(lo == hi, "%s: ranks disagreed (%s..%s)", what,
               h5c_status_string((h5c_status_t)lo),
               h5c_status_string((h5c_status_t)hi));
}

/*
 * Only ONE rank passes a bad index. Every rank must still come back with
 * H5C_ERR_INVALID_ARG, and none may enter a collective HDF5 call: if the
 * agreement were missing this would hang and the watchdog would fire.
 */
static void test_bad_connectivity(h5c_viz_t *viz)
{
    h5c_viz_mesh_t mesh;
    h5c_status_t   st;
    layout_t       lay;
    double        *nodes;
    int32_t       *conn;
    size_t         i, c, k;
    int            culprit = g_nprocs - 1;

    layout_of((size_t)(4 + g_me), (size_t)(2 + g_me), &lay);

    memset(&mesh, 0, sizeof mesh);
    mesh.kind             = H5C_VIZ_UNSTRUCTURED;
    mesh.name             = "rejected";
    mesh.topology         = "Triangle";
    mesh.nodes_per_element = NARROW_NPE;
    mesh.num_points       = lay.np;
    mesh.num_cells        = lay.nc;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    nodes = (double *)malloc(lay.np * 3 * sizeof *nodes);
    for (i = 0; i < lay.np; i++) {
        for (c = 0; c < 3; c++) {
            nodes[i * 3 + c] = node_c((int)c, g_me, i);
        }
    }
    H5C_CHECK(h5c_viz_write_nodes(viz, nodes, H5C_F64));
    free(nodes);

    conn = (int32_t *)malloc(lay.nc * NARROW_NPE * sizeof *conn);
    for (c = 0; c < lay.nc; c++) {
        for (k = 0; k < NARROW_NPE; k++) {
            conn[c * NARROW_NPE + k] = (int32_t)conn_local(c, k, lay.np);
        }
    }
    /* num_points itself is out of range: the classic 1-origin mistake. */
    if (g_me == culprit) {
        conn[0] = (int32_t)lay.np;
    }
    st = h5c_viz_write_connectivity(viz, conn, H5C_I32);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "out-of-range connectivity");

    /* A negative index is equally rejected, this time on rank 0. */
    if (g_me == culprit) {
        conn[0] = (int32_t)conn_local(0, 0, lay.np);
    }
    if (g_me == 0) {
        conn[0] = -1;
    }
    st = h5c_viz_write_connectivity(viz, conn, H5C_I32);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "negative connectivity");

    /* Valid input on the same mesh must still go through afterwards. */
    if (g_me == 0) {
        conn[0] = (int32_t)conn_local(0, 0, lay.np);
    }
    st = h5c_viz_write_connectivity(viz, conn, H5C_I32);
    H5C_CHECK(st);
    check_agreed(st, "valid connectivity after a rejection");
    free(conn);
}

/* Argument errors that need no mesh-wide data. */
static void test_bad_args(h5c_viz_t *viz)
{
    h5c_viz_mesh_t mesh;
    h5c_status_t   st;
    double         one[3] = { 1.0, 2.0, 3.0 };
    const void    *comps[3];

    comps[0] = one;
    comps[1] = one;
    comps[2] = one;

    memset(&mesh, 0, sizeof mesh);
    mesh.kind       = H5C_VIZ_POLYDATA;
    mesh.name       = "badargs";
    mesh.num_points = 1;
    mesh.num_cells  = 3;   /* a point cloud owns no cells */
    st = h5c_viz_begin_mesh(viz, &mesh);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "POLYDATA with cells");

    memset(&mesh, 0, sizeof mesh);
    mesh.kind = (h5c_viz_kind_t)77;
    st = h5c_viz_begin_mesh(viz, &mesh);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "unknown mesh kind");

    /* The failed begin_mesh must have left no mesh current. */
    st = h5c_viz_write_nodes(viz, one, H5C_F64);
    H5C_CHECK_FAILS(st, H5C_ERR_STATE);
    check_agreed(st, "write with no current mesh");

    memset(&mesh, 0, sizeof mesh);
    mesh.kind       = H5C_VIZ_POLYDATA;
    mesh.name       = "badargs";
    mesh.num_points = 1;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    /* A point cloud has neither connectivity nor cell data. */
    st = h5c_viz_write_connectivity(viz, one, H5C_I32);
    H5C_CHECK_FAILS(st, H5C_ERR_STATE);
    check_agreed(st, "connectivity on a point cloud");

    st = h5c_viz_write_cell_data(viz, "Nope", one, H5C_I32, 1);
    H5C_CHECK_FAILS(st, H5C_ERR_STATE);
    check_agreed(st, "cell data on a point cloud");

    st = h5c_viz_write_nodes(viz, one, H5C_I32);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "integer nodes");

    st = h5c_viz_write_point_data(viz, "Zero", one, H5C_F64, 0);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "ncomp 0");

    st = h5c_viz_write_point_data_comps(viz, "NullComps", NULL, H5C_F64, 3);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "NULL comps");

    st = h5c_viz_write_point_data(viz, "", one, H5C_F64, 1);
    H5C_CHECK_FAILS(st, H5C_ERR_INVALID_ARG);
    check_agreed(st, "empty field name");

    /* An unnamed component count is written, just without attribute_type. */
    H5C_CHECK(h5c_viz_write_point_data_comps(viz, "Odd", comps, H5C_F64, 2));

    /* Writing the same name twice is a mistake, not an overwrite. */
    H5C_CHECK(h5c_viz_write_nodes(viz, one, H5C_F64));
    st = h5c_viz_write_nodes(viz, one, H5C_F64);
    H5C_CHECK_FAILS(st, H5C_ERR_EXISTS);
    check_agreed(st, "nodes written twice");
}

static void check_bad_args(h5c_file_t *f)
{
    H5C_ASSERT(h5c_exists(f, "/badargs/point_data/Odd"),
               "the 2-component field should have been written");
    H5C_ASSERT(!h5c_attr_exists(f, "/badargs/point_data/Odd",
                                "attribute_type"),
               "a 2-component field has no XDMF attribute_type");
    H5C_ASSERT(!h5c_exists(f, "/badargs/point_data/Zero"),
               "the ncomp 0 field must not exist");
    H5C_ASSERT(h5c_exists(f, "/rejected/geometry/connectivity"),
               "the connectivity retried with valid indices is missing");
    H5C_ASSERT(!h5c_exists(f, "/badargs/cell_data"),
               "a point cloud must have no cell_data group");
}

/* ------------------------------------------------------------------ */
/* attribute_type mapping                                              */
/* ------------------------------------------------------------------ */

static void test_attribute_type(void)
{
    const char *got;

    H5C_ASSERT(strcmp(h5c_viz_attribute_type(1), "Scalar") == 0,
               "ncomp 1 is not Scalar");
    H5C_ASSERT(strcmp(h5c_viz_attribute_type(3), "Vector") == 0,
               "ncomp 3 is not Vector");
    H5C_ASSERT(strcmp(h5c_viz_attribute_type(6), "Tensor6") == 0,
               "ncomp 6 is not Tensor6");
    H5C_ASSERT(strcmp(h5c_viz_attribute_type(9), "Tensor") == 0,
               "ncomp 9 is not Tensor");
    got = h5c_viz_attribute_type(0);
    H5C_ASSERT(got == NULL, "ncomp 0 mapped to '%s'", got != NULL ? got : "");
    got = h5c_viz_attribute_type(2);
    H5C_ASSERT(got == NULL, "ncomp 2 mapped to '%s'", got != NULL ? got : "");
    got = h5c_viz_attribute_type(4);
    H5C_ASSERT(got == NULL, "ncomp 4 mapped to '%s'", got != NULL ? got : "");
}

/* ------------------------------------------------------------------ */

static void check_root(h5c_file_t *f)
{
    int32_t scheme = -1;
    double  time = -1.0;

    H5C_CHECK(h5c_read_attr_scalar(f, "/", "scheme_version", &scheme,
                                   H5C_I32));
    H5C_ASSERT(scheme == H5C_SCHEME_VERSION,
               "scheme_version: got %d, want %d", (int)scheme,
               (int)H5C_SCHEME_VERSION);
    H5C_CHECK(h5c_read_attr_scalar(f, "/", "time", &time, H5C_F64));
    H5C_ASSERT(time == TIME_VALUE, "time: got %.17g, want %.17g", time,
               (double)TIME_VALUE);
}

int main(int argc, char **argv)
{
    h5c_viz_t  *viz = NULL;
    h5c_file_t *f = NULL;
    int         total = 0;
    int         i;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_me);
    MPI_Comm_size(MPI_COMM_WORLD, &g_nprocs);

    signal(SIGALRM, on_watchdog);
    alarm(WATCHDOG_SECONDS);

    H5C_CHECK(h5c_init());
    test_attribute_type();

    H5C_CHECK(h5c_viz_open(PATH, TIME_VALUE, MPI_COMM_WORLD, MPI_INFO_NULL,
                           &viz));
    if (viz != NULL) {
        write_fluid(viz);
        write_dust(viz);
        write_sparse(viz);
        for (i = 0; i < 4; i++) {
            write_narrow(viz, &g_narrow[i]);
        }
        /* Everything above must have gone through cleanly. */
        H5C_CHECK(h5c_viz_status(viz));

        /* From here on failures are expected, so the sticky status is used. */
        test_bad_connectivity(viz);
        test_bad_args(viz);
        H5C_CHECK_FAILS(h5c_viz_status(viz), H5C_ERR_INVALID_ARG);

        H5C_CHECK(h5c_viz_close(viz));

        /*
         * Read back through the parallel entry point: these are plain
         * datasets, not the distributed __partition__ layout, so every rank
         * reads the whole dataset independently and checks its own slice.
         */
        H5C_CHECK(h5c_popen(PATH, H5C_READ, &f));
        if (f != NULL) {
            check_root(f);
            check_fluid(f);
            check_dust(f);
            check_sparse(f);
            for (i = 0; i < 4; i++) {
                check_narrow(f, &g_narrow[i]);
            }
            check_bad_args(f);
            H5C_CHECK(h5c_close(f));
        }
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
    return H5C_TEST_SUMMARY("test_pviz");
}
