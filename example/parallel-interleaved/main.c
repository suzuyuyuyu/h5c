/*
 * Distributed vector field: the combination a solver actually needs.
 *
 * Components live in separate arrays (u, v, w) because that is what suits a
 * GPU kernel, the points are split across ranks, and the file has to come out
 * as one [total_n, 3] dataset so that XDMF sees a vector.
 *
 * h5c does all three at once: it packs each rank's components into a
 * contiguous staging buffer, writes them at that rank's offset, and records
 * the rank boundaries beside the data.
 *
 * MUST be run through the batch system, never on a login node:
 *     sbatch example/run-parallel-example.sh
 */
#include <h5c/h5c_mpi.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NCOMP 3

int main(int argc, char **argv)
{
    const char *path = "example_parallel_interleaved.h5";

    int    me = 0, nprocs = 1;
    size_t nlocal, offset = 0, mine = 0;
    double *u = NULL, *v = NULL, *w = NULL;
    double *gu = NULL, *gv = NULL, *gw = NULL;
    const double *comps[NCOMP];
    double       *out[NCOMP];
    h5c_file_t   *f = NULL;
    h5c_status_t  st;
    size_t        i;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /*
     * Rank 0 deliberately owns nothing when there is more than one rank.
     * An empty subdomain is legal: the rank contributes no rows but still
     * takes part in every collective call, and its component pointers may
     * be NULL.
     */
    nlocal = (nprocs > 1 && me == 0) ? 0u : (size_t)(2 + me);

    if (nlocal > 0) {
        u = (double *)malloc(nlocal * sizeof *u);
        v = (double *)malloc(nlocal * sizeof *v);
        w = (double *)malloc(nlocal * sizeof *w);
        gu = (double *)malloc(nlocal * sizeof *gu);
        gv = (double *)malloc(nlocal * sizeof *gv);
        gw = (double *)malloc(nlocal * sizeof *gw);
        for (i = 0; i < nlocal; i++) {
            /* Every value encodes its rank and position, so a mis-offset
               write cannot look right by accident. */
            u[i] = 1000.0 * me + (double)i;
            v[i] = 2000.0 * me + (double)i;
            w[i] = 3000.0 * me + (double)i;
        }
    }
    comps[0] = u; comps[1] = v; comps[2] = w;
    out[0] = gu;  out[1] = gv;  out[2] = gw;

    /* ---- write ------------------------------------------------------ */

    if ((st = h5c_popen(path, H5C_TRUNCATE, &f)) != H5C_OK) {
        if (me == 0) {
            fprintf(stderr, "popen: %s\n", h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    st = h5c_pwrite_interleaved(f, "/fields/velocity",
                                (const void *const *)comps,
                                NCOMP, nlocal, H5C_F64, H5C_WRITE_DEFAULT);
    if (st != H5C_OK) {
        if (me == 0) {
            fprintf(stderr, "pwrite_interleaved: %s (%s)\n",
                    h5c_status_string(st), h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /*
     * The dataset lives at "<path>/data"; the group also holds __partition__.
     * XDMF needs this attribute to treat the three components as a vector.
     */
    h5c_write_attr_str(f, "/fields/velocity/data", "attribute_type", "Vector");

    h5c_poffset(f, "/fields/velocity", &offset, &mine);
    printf("rank %d: %lu points at rows [%lu, %lu)\n", me,
           (unsigned long)mine, (unsigned long)offset,
           (unsigned long)(offset + mine));
    fflush(stdout);

    h5c_close(f);

    /* ---- read back -------------------------------------------------- */

    if ((st = h5c_popen(path, H5C_READ, &f)) != H5C_OK) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    st = h5c_pread_interleaved(f, "/fields/velocity", (void *const *)out,
                               NCOMP, nlocal, H5C_F64);
    if (st != H5C_OK) {
        if (me == 0) {
            fprintf(stderr, "pread_interleaved: %s (%s)\n",
                    h5c_status_string(st), h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (i = 0; i < nlocal; i++) {
        if (gu[i] != u[i] || gv[i] != v[i] || gw[i] != w[i]) {
            fprintf(stderr, "rank %d: point %lu differs\n",
                    me, (unsigned long)i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Rank 0 shows that the file is one interleaved [total_n, 3] dataset,
       not three separate blocks. Every rank must reach the collective calls
       below, so only the PRINTING is guarded. */
    {
        h5c_dataset_info_t global;
        h5c_pdataset_info(f, "/fields/velocity", NULL, &global);

        if (me == 0) {
            double *flat = (double *)malloc(global.count * sizeof *flat);
            h5c_file_t *serial = NULL;

            printf("global shape: {%lu, %lu} from %d ranks\n",
                   (unsigned long)global.dims[0],
                   (unsigned long)global.dims[1], nprocs);

            /* Reading the whole thing is a plain, undistributed read. */
            if (h5c_open(path, H5C_READ, &serial) == H5C_OK) {
                h5c_read_f64(serial, "/fields/velocity/data", flat,
                             2, global.dims);
                printf("as stored:");
                for (i = 0; i < global.count && i < 12; i++) {
                    printf(" %g", flat[i]);
                }
                printf("%s\n", (global.count > 12) ? " ..." : "");
                h5c_close(serial);
            }
            free(flat);
        }
    }

    h5c_close(f);

    free(u); free(v); free(w);
    free(gu); free(gv); free(gw);

    if (me == 0) {
        printf("wrote and read %s\n", path);
    }
    MPI_Finalize();
    return 0;
}
