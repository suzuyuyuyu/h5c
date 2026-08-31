/*
 * Distributed read and write.
 *
 * Each rank owns a slice of the slowest-varying axis. h5c concatenates the
 * slices into one dataset and records the rank boundaries beside it, so the
 * result is a single file that any number of ranks can read back.
 *
 * MUST be run through the batch system, never on a login node:
 *     sbatch scripts/run-mpi-tests.sh    (or your own job script)
 */
#include <h5c/h5c_mpi.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    const char *path = "example_parallel.h5";

    int    me = 0, nprocs = 1;
    size_t nlocal, dims[2], offset = 0, mine = 0;
    double *local = NULL, *got = NULL;
    h5c_file_t *f = NULL;
    h5c_status_t st;
    size_t i;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Deliberately unequal: rank r owns 2 + r rows. A rank may own 0. */
    nlocal  = (size_t)(2 + me);
    dims[0] = nlocal;   /* the split axis */
    dims[1] = 3;        /* must be the same on every rank */

    local = (double *)malloc(nlocal * 3 * sizeof *local);
    got   = (double *)malloc(nlocal * 3 * sizeof *got);
    for (i = 0; i < nlocal * 3; i++) {
        local[i] = 100.0 * me + (double)i;
    }

    /* Every call below is collective: all ranks, same order, same path. */
    if ((st = h5c_popen(path, H5C_TRUNCATE, &f)) != H5C_OK) {
        if (me == 0) {
            fprintf(stderr, "popen: %s\n", h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    st = h5c_pwrite(f, "/coords", local, H5C_F64, 2, dims, H5C_WRITE_DEFAULT);
    if (st != H5C_OK) {
        /* Argument checks are agreed across ranks before any collective HDF5
           call, so every rank lands here together rather than deadlocking. */
        if (me == 0) {
            fprintf(stderr, "pwrite: %s (%s)\n", h5c_status_string(st),
                    h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Where did this rank's block end up? No need to open __partition__. */
    h5c_poffset(f, "/coords", &offset, &mine);
    printf("rank %d: rows [%lu, %lu)\n", me,
           (unsigned long)offset, (unsigned long)(offset + mine));
    fflush(stdout);

    h5c_close(f);

    /* ---- read back ------------------------------------------------ */

    if ((st = h5c_popen(path, H5C_READ, &f)) != H5C_OK) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* The stored boundaries must match the current rank count; h5c checks
       that before moving any data. */
    if ((st = h5c_pread(f, "/coords", got, H5C_F64, 2, dims)) != H5C_OK) {
        if (me == 0) {
            fprintf(stderr, "pread: %s (%s)\n", h5c_status_string(st),
                    h5c_last_error()->message);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (i = 0; i < nlocal * 3; i++) {
        if (got[i] != local[i]) {
            fprintf(stderr, "rank %d: element %lu differs\n",
                    me, (unsigned long)i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    {
        h5c_dataset_info_t global;
        h5c_pdataset_info(f, "/coords", NULL, &global);
        if (me == 0) {
            printf("global shape: {%lu, %lu} from %d ranks\n",
                   (unsigned long)global.dims[0],
                   (unsigned long)global.dims[1], nprocs);
        }
    }

    h5c_close(f);
    free(local);
    free(got);

    if (me == 0) {
        printf("wrote and read %s\n", path);
    }
    MPI_Finalize();
    return 0;
}
