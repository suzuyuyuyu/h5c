#include "h5c_internal.h"

#include <stdarg.h>
#include <string.h>

static H5C_THREAD_LOCAL h5c_error_t g_last = { H5C_OK, 0, { '\0' } };

static int g_verbosity = 0;
static int g_initialised = 0;

const char *h5c_status_string(h5c_status_t status)
{
    switch (status) {
    case H5C_OK:                 return "ok";
    case H5C_ERR_INVALID_ARG:    return "invalid argument";
    case H5C_ERR_NOT_FOUND:      return "not found";
    case H5C_ERR_SHAPE_MISMATCH: return "shape mismatch";
    case H5C_ERR_TYPE_MISMATCH:  return "type mismatch";
    case H5C_ERR_EXISTS:         return "already exists";
    case H5C_ERR_HDF5:           return "HDF5 error";
    case H5C_ERR_MPI:            return "MPI error";
    case H5C_ERR_NOMEM:          return "out of memory";
    case H5C_ERR_STATE:          return "invalid state";
    case H5C_ERR_UNSUPPORTED:    return "unsupported";
    default:                     return "unknown status";
    }
}

const h5c_error_t *h5c_last_error(void)
{
    return &g_last;
}

static void store(h5c_status_t status, long herr, const char *fmt, va_list ap)
{
    g_last.status   = status;
    g_last.hdf5_err = herr;
    if (fmt != NULL) {
        vsnprintf(g_last.message, sizeof g_last.message, fmt, ap);
    } else {
        g_last.message[0] = '\0';
    }
}

h5c_status_t h5c__fail(h5c_status_t status, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    store(status, 0, fmt, ap);
    va_end(ap);
    return status;
}

h5c_status_t h5c__fail_hdf5(long herr, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    store(H5C_ERR_HDF5, herr, fmt, ap);
    va_end(ap);
    return H5C_ERR_HDF5;
}

h5c_status_t h5c__record(h5c_file_t *file, h5c_status_t status)
{
    if (file != NULL && file->sticky == H5C_OK && status != H5C_OK) {
        file->sticky = status;
    }
    return status;
}

void h5c_set_error_verbosity(int level)
{
    g_verbosity = level;
    if (g_initialised) {
        if (level > 0) {
            H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr);
        } else {
            H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
        }
    }
}

h5c_status_t h5c__ensure_init(void)
{
    if (g_initialised) {
        return H5C_OK;
    }
    if (H5open() < 0) {
        return h5c__fail_hdf5(-1, "H5open failed");
    }
    g_initialised = 1;
    /* Apply the configured verbosity now that HDF5 is up. */
    h5c_set_error_verbosity(g_verbosity);
    return H5C_OK;
}

h5c_status_t h5c_init(void)
{
    return h5c__ensure_init();
}

void h5c_finalize(void)
{
    if (!g_initialised) {
        return;
    }
    h5c__type_cleanup();
    H5close();
    g_initialised = 0;
}
