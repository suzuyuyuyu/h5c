/*
 * A time series a post-processor can actually consume.
 *
 * Writes one HDF5 file per step, each holding two meshes: a tetrahedral grid
 * with a pressure wave and a swirling velocity, and a particle cloud that
 * settles and spreads. That is the shape of a real solver's visualization
 * output, and it is what h5xdmf turns into XDMF3 for ParaView.
 *
 * The point of this example is that the files h5c writes are the SAME layout
 * h5fortran writes, so the same Python tooling reads either.
 *
 * MUST be run through the batch system, never on a login node:
 *     sbatch example/visualization/run.sh
 */
#include <h5c/h5c_viz.h>

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NSTEPS   5
#define NCELLS_X 4          /* tetra cells this rank contributes per step */
#define NPARTS   6          /* particles per rank */

static int g_me = 0, g_nprocs = 1;

/* Fails loudly and identically on every rank: these calls are collective. */
static void must(h5c_status_t st, const char *what)
{
    if (st != H5C_OK) {
        if (g_me == 0) {
            fprintf(stderr, "%s: %s (%s)\n", what, h5c_status_string(st),
                    h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

/*
 * One tetrahedron per cell, each with its own four nodes. Duplicating the
 * shared nodes keeps the example short; a real solver would share them within
 * a rank and let only the rank boundaries duplicate.
 */
static void build_tetra(size_t ncells, double t,
                        double **nodes, int32_t **conn,
                        double **pressure, double **velocity,
                        int32_t **subdomain)
{
    const size_t npoints = ncells * 4;
    size_t c, k;

    *nodes     = (double *)malloc(npoints * 3 * sizeof(double));
    *conn      = (int32_t *)malloc(ncells * 4 * sizeof(int32_t));
    *pressure  = (double *)malloc(npoints * sizeof(double));
    *velocity  = (double *)malloc(npoints * 3 * sizeof(double));
    *subdomain = (int32_t *)malloc(ncells * sizeof(int32_t));

    for (c = 0; c < ncells; c++) {
        /* Offset each cell so the mesh spreads along x, ranks along y. */
        const double ox = (double)c;
        const double oy = (double)g_me;

        static const double unit[4][3] = {
            { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }
        };
        for (k = 0; k < 4; k++) {
            const size_t p = c * 4 + k;
            const double x = ox + unit[k][0];
            const double y = oy + unit[k][1];
            const double z = unit[k][2];

            (*nodes)[p * 3 + 0] = x;
            (*nodes)[p * 3 + 1] = y;
            (*nodes)[p * 3 + 2] = z;

            /* A wave travelling in x, so the series actually animates. */
            (*pressure)[p] = sin(x - 2.0 * t) * cos(0.5 * y);
            /* Swirl about the z axis. */
            (*velocity)[p * 3 + 0] = -y;
            (*velocity)[p * 3 + 1] =  x;
            (*velocity)[p * 3 + 2] =  0.1 * sin(t);

            /* Connectivity is RANK-LOCAL and 0-origin; h5c adds the offset. */
            (*conn)[c * 4 + k] = (int32_t)p;
        }
        (*subdomain)[c] = g_me;
    }
}

static void write_step(int step, double t)
{
    char path[256];
    h5c_viz_t *viz = NULL;
    h5c_viz_mesh_t mesh;

    /* Zero-padded so lexical order is time order. */
    snprintf(path, sizeof path, "result/seq%06d.h5", step);

    must(h5c_viz_open(path, t, MPI_COMM_WORLD, MPI_INFO_NULL, &viz),
         "viz_open");

    /* ---- the tetrahedral grid ------------------------------------- */
    {
        const size_t ncells  = NCELLS_X;
        const size_t npoints = ncells * 4;
        double *nodes, *pressure, *velocity;
        int32_t *conn, *subdomain;

        build_tetra(ncells, t, &nodes, &conn, &pressure, &velocity, &subdomain);

        memset(&mesh, 0, sizeof mesh);
        mesh.kind              = H5C_VIZ_UNSTRUCTURED;
        mesh.name              = "fluid";
        mesh.topology          = "Tetrahedron";
        mesh.nodes_per_element = 4;
        mesh.num_points        = npoints;
        mesh.num_cells         = ncells;
        must(h5c_viz_begin_mesh(viz, &mesh), "begin_mesh fluid");

        must(h5c_viz_write_nodes(viz, nodes, H5C_F64), "write_nodes");
        must(h5c_viz_write_connectivity(viz, conn, H5C_I32), "write_conn");

        must(h5c_viz_write_point_data(viz, "Pressure", pressure, H5C_F64, 1),
             "point Pressure");
        must(h5c_viz_write_point_data(viz, "Velocity", velocity, H5C_F64, 3),
             "point Velocity");
        must(h5c_viz_write_cell_data(viz, "SubdomainID", subdomain, H5C_I32, 1),
             "cell SubdomainID");

        free(nodes); free(conn); free(pressure); free(velocity);
        free(subdomain);
    }

    /* ---- the particle cloud --------------------------------------- */
    {
        double x[NPARTS], y[NPARTS], z[NPARTS], radius[NPARTS];
        const void *xyz[3];
        size_t i;

        for (i = 0; i < NPARTS; i++) {
            const double s = (double)i / (double)NPARTS;
            x[i] = 0.5 + 3.0 * s + 0.3 * sin(t + s);
            y[i] = (double)g_me + 0.5 + 0.2 * cos(t + s);
            z[i] = 2.0 - 0.15 * t * (1.0 + s);   /* settling */
            radius[i] = 0.05 + 0.02 * s;
        }
        xyz[0] = x; xyz[1] = y; xyz[2] = z;

        memset(&mesh, 0, sizeof mesh);
        mesh.kind       = H5C_VIZ_POLYDATA;
        mesh.name       = "particles";
        mesh.num_points = NPARTS;
        mesh.num_cells  = 0;
        must(h5c_viz_begin_mesh(viz, &mesh), "begin_mesh particles");

        /* Coordinates held as separate arrays: no packing by the caller. */
        must(h5c_viz_write_nodes_comps(viz, xyz, H5C_F64), "write_nodes_comps");
        must(h5c_viz_write_point_data(viz, "Radius", radius, H5C_F64, 1),
             "point Radius");
    }

    must(h5c_viz_status(viz), "viz writes");
    must(h5c_viz_close(viz), "viz_close");

    if (g_me == 0) {
        printf("wrote %s (t = %.3f)\n", path, t);
        fflush(stdout);
    }
}

int main(int argc, char **argv)
{
    int step;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_me);
    MPI_Comm_size(MPI_COMM_WORLD, &g_nprocs);

    for (step = 0; step < NSTEPS; step++) {
        write_step(step, 0.25 * (double)step);
    }

    if (g_me == 0) {
        printf("\n%d steps from %d ranks. Now generate XDMF:\n"
               "  cd ../h5fortran/postprocess && uv sync\n"
               "  uv run h5xdmf \"<this dir>/result/seq*.h5\" \\\n"
               "      --metadata <this dir>/result/metadata.h5 \\\n"
               "      --outdir <this dir>/result\n"
               "then open result/fluid.xdmf and result/particles.xdmf"
               " in ParaView.\n", NSTEPS, g_nprocs);
    }
    MPI_Finalize();
    return 0;
}
