/* Parallel visualization opening; all subsequent operations are shared. */
#ifndef H5C_VIZ_MPI_H
#define H5C_VIZ_MPI_H

#include <mpi.h>

#include "h5c/h5c_viz.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Collective over comm, which must remain valid until h5c_viz_close(). */
h5c_status_t h5c_viz_popen(const char* path, double time, MPI_Comm comm, MPI_Info info, h5c_viz_t** out);

#ifdef __cplusplus
}
#endif
#endif
