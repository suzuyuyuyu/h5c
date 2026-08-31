/*
 * Serial read and write.
 *
 * Covers the shape of ordinary h5c code: open, write a few datasets, attach
 * attributes, then read back by asking for the shape first.
 */
#include <h5c/h5c.h>

#include <stdio.h>
#include <stdlib.h>

/* Reports the recorded message, which says more than the status alone. */
static int failed(const char *what, h5c_status_t st)
{
    if (st == H5C_OK) {
        return 0;
    }
    fprintf(stderr, "%s: %s (%s)\n", what, h5c_status_string(st),
            h5c_last_error()->message);
    return 1;
}

int main(void)
{
    const char *path = "example_serial.h5";

    /* dims are ROW-MAJOR: the last one varies fastest.
       This is the C view of what Fortran declares as a(2, 3). */
    const double values[6] = { 1, 2, 3, 4, 5, 6 };
    const size_t dims[2]   = { 3, 2 };

    const h5c_bool_t flags[4] = { H5C_TRUE, H5C_FALSE, H5C_FALSE, H5C_TRUE };
    const size_t     fdims[2] = { 2, 2 };

    h5c_file_t         *f = NULL;
    h5c_dataset_info_t  info;
    double             *got = NULL;
    char               *units = NULL;
    size_t              i;

    /* ---- write ---------------------------------------------------- */

    if (failed("open", h5c_open(path, H5C_TRUNCATE, &f))) {
        return 1;
    }

    /* Intermediate groups are created as needed. */
    h5c_write_f64(f, "/mesh/coords", values, 2, dims);
    h5c_write_f64_scalar(f, "/time", 0.125);
    h5c_write_bool(f, "/mesh/active", flags, 2, fdims);
    h5c_write_string(f, "/title", "example field", H5C_WRITE_DEFAULT);

    /* Attributes are written separately from the data they annotate. */
    h5c_write_attr_str(f, "/mesh/coords", "units", "m");
    h5c_write_attr_str(f, "/", "created_by", "h5c example");

    /* One check at the end is enough: the file keeps the FIRST failure. */
    if (failed("write", h5c_file_status(f))) {
        h5c_close(f);
        return 1;
    }
    /* h5c_close reports only whether closing worked, so a failure here
       really does mean the file may not have been flushed. */
    if (failed("close", h5c_close(f))) {
        return 1;
    }

    /* ---- read back ------------------------------------------------ */

    if (failed("reopen", h5c_open(path, H5C_READ, &f))) {
        return 1;
    }

    /* Ask for the shape, then allocate: h5c never allocates behind your back
       unless you call the *_alloc form. */
    if (failed("info", h5c_dataset_info(f, "/mesh/coords", &info))) {
        h5c_close(f);
        return 1;
    }
    printf("/mesh/coords: rank=%d dims={%lu, %lu} count=%lu\n",
           info.rank, (unsigned long)info.dims[0],
           (unsigned long)info.dims[1], (unsigned long)info.count);

    got = (double *)malloc(info.count * sizeof *got);
    if (got == NULL) {
        h5c_close(f);
        return 1;
    }
    if (failed("read", h5c_read_f64(f, "/mesh/coords", got,
                                    info.rank, info.dims))) {
        free(got);
        h5c_close(f);
        return 1;
    }
    printf("values:");
    for (i = 0; i < info.count; i++) {
        printf(" %g", got[i]);
    }
    putchar('\n');
    free(got);

    /* Strings are allocated by h5c and released with h5c_free_string(). */
    if (h5c_read_attr_str(f, "/mesh/coords", "units", &units) == H5C_OK) {
        printf("units: %s\n", units);
        h5c_free_string(units);
    }

    /* A missing dataset is an ordinary status, not a crash. */
    {
        double unused = 0.0;
        h5c_status_t st = h5c_read_f64_scalar(f, "/not/here", &unused);
        printf("missing dataset -> %s\n", h5c_status_string(st));
        h5c_file_clear_status(f);   /* we handled it; drop it from the file */
    }

    h5c_close(f);
    printf("wrote and read %s\n", path);
    return 0;
}
