/*
 * Vector fields held as separate component arrays.
 *
 * Solvers usually keep u, v, w apart, but XDMF wants one [n, 3] dataset so
 * that ParaView sees a vector (velocity magnitude, tensor invariants). h5c
 * does the interleaving, so the solver keeps its own layout.
 */
#include <h5c/h5c.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4

int main(void) {
    const char* path = "example_serial_interleaved.h5";

    /* Separate components, as a solver would hold them. */
    double u[N] = {1, 2, 3, 4};
    double v[N] = {10, 20, 30, 40};
    double w[N] = {100, 200, 300, 400};

    const double* comps[3];
    double* out[3];
    double gu[N], gv[N], gw[N];
    double flat[N * 3];
    double only_v[N];

    h5c_file_t* f = NULL;
    h5c_status_t st;
    size_t i;

    comps[0] = u;
    comps[1] = v;
    comps[2] = w;

    if ((st = h5c_open(path, H5C_TRUNCATE, &f)) != H5C_OK) {
        fprintf(stderr, "open: %s\n", h5c_last_error()->message);
        return 1;
    }

    /* One call; h5c packs the components into a contiguous [N, 3] buffer.
       Writing packs rather than striding because a strided collective write
       forces read-modify-write in MPI-IO. */
    h5c_write_interleaved_f64(f, "/fields/velocity", comps, 3, N);

    /* Vector fields need this attribute for the XDMF post-processing. */
    h5c_write_attr_str(f, "/fields/velocity", "attribute_type", "Vector");

    if ((st = h5c_file_status(f)) != H5C_OK) {
        fprintf(stderr, "write: %s (%s)\n", h5c_status_string(st), h5c_last_error()->message);
        h5c_close(f);
        return 1;
    }
    h5c_close(f);

    /* ---- read back ------------------------------------------------ */

    if ((st = h5c_open(path, H5C_READ, &f)) != H5C_OK) {
        fprintf(stderr, "reopen: %s\n", h5c_last_error()->message);
        return 1;
    }

    /* Scattered straight back into separate arrays. */
    out[0] = gu;
    out[1] = gv;
    out[2] = gw;
    h5c_read_interleaved_f64(f, "/fields/velocity", out, 3, N);
    printf("components:\n");
    for (i = 0; i < N; i++) {
        printf("  %g %g %g\n", gu[i], gv[i], gw[i]);
    }

    /* The file really is interleaved: reading it as a plain [N, 3] dataset
       gives u0 v0 w0 u1 v1 w1 ... */
    {
        const size_t dims[2] = {N, 3};
        h5c_read_f64(f, "/fields/velocity", flat, 2, dims);
        printf("as stored:");
        for (i = 0; i < N * 3; i++) {
            printf(" %g", flat[i]);
        }
        putchar('\n');
    }

    /* Only one component? Then only its bytes move: a strided read has no
       read-modify-write penalty, unlike a strided write. */
    h5c_read_component_f64(f, "/fields/velocity", only_v, 1, N);
    printf("v only:");
    for (i = 0; i < N; i++) {
        printf(" %g", only_v[i]);
    }
    putchar('\n');

    h5c_close(f);
    return 0;
}
