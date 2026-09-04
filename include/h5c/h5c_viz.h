/*
 * h5c — visualization HDF5 writer (scheme_version 1).
 *
 * Writes the layout the h5xdmf post-processor reads, identical to the one
 * h5fortran's t_phdf5_writer produces, so the same Python tooling turns either
 * into XDMF3. No XML is generated here.
 *
 * ---------------------------------------------------------------------------
 * File layout
 * ---------------------------------------------------------------------------
 *
 *   /                              attrs: scheme_version=1, time=<f64>
 *   /<mesh>/                       attrs: topology_type, nodes_per_element
 *   /<mesh>/geometry/nodes         (total_points, 3)
 *   /<mesh>/geometry/connectivity  (total_cells, nodes_per_element)
 *   /<mesh>/point_data/<field>     (total_points[, ncomp])
 *   /<mesh>/cell_data/<field>      (total_cells[, ncomp])
 *
 * connectivity is absent for a point cloud (H5C_VIZ_POLYDATA). Fields with
 * ncomp > 1 also carry attribute_type = Vector (3), Tensor6 (6) or Tensor (9);
 * ncomp == 1 is a plain 1-D dataset and XDMF treats it as a Scalar.
 *
 * Several meshes may live in one file: call h5c_viz_begin_mesh() again with a
 * different name. h5xdmf emits one .xdmf per mesh group.
 *
 * ---------------------------------------------------------------------------
 * What each rank contributes
 * ---------------------------------------------------------------------------
 *
 * num_cells is the number of cells this rank OWNS; ghost cells must not be
 * written or they appear twice. num_points covers every local node those
 * cells reference, so a node shared across a rank boundary is written once per
 * rank and the duplicates are expected -- ParaView's "Clean to Grid" merges
 * them if that matters.
 *
 * Connectivity is given as RANK-LOCAL, 0-origin node indices into this rank's
 * own nodes. The writer adds this rank's node offset to produce file-global
 * ids, so the caller needs no communication to number its nodes.
 *
 * ---------------------------------------------------------------------------
 * Collective discipline
 * ---------------------------------------------------------------------------
 *
 * Every function here is collective over the file's communicator: all ranks
 * must call it, in the same order, with the same mesh and field names. Point
 * and cell counts may differ per rank and may be 0.
 */
#ifndef H5C_VIZ_H
#define H5C_VIZ_H

#include "h5c/h5c_mpi.h"
#include "h5c/h5c_version.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Defaults applied when a mesh leaves the corresponding field unset. */
#define H5C_VIZ_DEFAULT_UGRID_NAME    "ugrid"
#define H5C_VIZ_DEFAULT_POLYDATA_NAME "polydata"
#define H5C_VIZ_DEFAULT_TOPOLOGY      "Hexahedron"
#define H5C_VIZ_DEFAULT_NODES_PER_ELEM 8

typedef struct h5c_viz h5c_viz_t;

typedef enum h5c_viz_kind {
    /* Cells plus connectivity. */
    H5C_VIZ_UNSTRUCTURED = 1,
    /* A point cloud: nodes only, no connectivity, no cell data. */
    H5C_VIZ_POLYDATA     = 2
} h5c_viz_kind_t;

/*
 * Describes one mesh. Zero-initialise and fill in what matters; NULL and 0
 * take the documented defaults.
 *
 *     h5c_viz_mesh_t mesh = {0};
 *     mesh.kind       = H5C_VIZ_UNSTRUCTURED;
 *     mesh.name       = "fluid";
 *     mesh.topology   = "Tetrahedron";
 *     mesh.nodes_per_element = 4;
 *     mesh.num_points = my_points;
 *     mesh.num_cells  = my_cells;
 */
typedef struct h5c_viz_mesh {
    h5c_viz_kind_t kind;
    const char *name;              /* group under "/"; NULL -> kind default */
    const char *topology;          /* XDMF topology_type; NULL -> Hexahedron */
    int    nodes_per_element;      /* 0 -> 8; ignored for POLYDATA */
    size_t num_points;             /* THIS rank's node count; may be 0 */
    size_t num_cells;              /* THIS rank's owned cells; 0 for POLYDATA */
} h5c_viz_mesh_t;

/* ------------------------------------------------------------------ */
/* lifecycle                                                           */
/* ------------------------------------------------------------------ */

/*
 * Creates `path`, truncating it, and stamps the root attributes
 * scheme_version = H5C_SCHEME_VERSION and time = `time`.
 *
 * One file per time step is the intended usage; `time` is what places it in
 * the series. `comm` may be MPI_COMM_WORLD and `info` may be MPI_INFO_NULL.
 */
h5c_status_t h5c_viz_open(const char *path, double time,
                          MPI_Comm comm, MPI_Info info, h5c_viz_t **out);

/*
 * Closes the file and frees the handle, which is invalid afterwards either
 * way. The status describes only the close, as for h5c_close(); read
 * h5c_viz_status() first if you want to know whether anything failed earlier.
 */
h5c_status_t h5c_viz_close(h5c_viz_t *viz);

/* First non-OK status recorded on this writer, or H5C_OK. */
h5c_status_t h5c_viz_status(const h5c_viz_t *viz);

/* ------------------------------------------------------------------ */
/* meshes                                                              */
/* ------------------------------------------------------------------ */

/*
 * Starts (or re-selects) a mesh and makes it the target of the write calls
 * below. Creates the group, its geometry/point_data/cell_data subgroups and
 * the topology attributes, and gathers the per-rank point and cell counts so
 * that offsets are known.
 *
 * Calling it again with a different name adds another mesh to the same file.
 */
h5c_status_t h5c_viz_begin_mesh(h5c_viz_t *viz, const h5c_viz_mesh_t *mesh);

/* This rank's offsets within the current mesh, once begin_mesh has run. */
h5c_status_t h5c_viz_offsets(const h5c_viz_t *viz,
                             size_t *point_offset, size_t *cell_offset);

/* ------------------------------------------------------------------ */
/* geometry                                                            */
/* ------------------------------------------------------------------ */

/*
 * Node coordinates, written to <mesh>/geometry/nodes as (total_points, 3).
 *
 * The first form takes one interleaved buffer of num_points * 3 values
 * (x0 y0 z0 x1 y1 z1 ...); the second takes three separate coordinate arrays,
 * which is how a solver that keeps x, y and z apart already holds them.
 * `type` must be H5C_F32 or H5C_F64.
 */
h5c_status_t h5c_viz_write_nodes(h5c_viz_t *viz, const void *nodes,
                                 h5c_type_t type);

h5c_status_t h5c_viz_write_nodes_comps(h5c_viz_t *viz,
                                       const void *const *xyz,
                                       h5c_type_t type);

/*
 * Cell connectivity, written to <mesh>/geometry/connectivity as
 * (total_cells, nodes_per_element).
 *
 * `conn` holds num_cells * nodes_per_element RANK-LOCAL 0-origin node
 * indices; the writer adds this rank's node offset. Any index outside
 * [0, num_points) is rejected with H5C_ERR_INVALID_ARG, which catches the
 * common mistake of passing 1-origin or already-global ids.
 *
 * `type` may be H5C_I8, H5C_I16, H5C_I32 or H5C_I64. Not valid for POLYDATA.
 */
h5c_status_t h5c_viz_write_connectivity(h5c_viz_t *viz, const void *conn,
                                        h5c_type_t type);

/* ------------------------------------------------------------------ */
/* fields                                                              */
/* ------------------------------------------------------------------ */

/*
 * Point and cell fields, written to <mesh>/point_data/<name> and
 * <mesh>/cell_data/<name>.
 *
 * ncomp == 1 gives a 1-D dataset; ncomp of 3, 6 or 9 gives (n, ncomp) plus
 * the matching attribute_type. Other values are accepted as (n, ncomp) but
 * carry no attribute_type, since XDMF has no name for them.
 *
 * For ncomp == 6 the XDMF component order is XX, XY, XZ, YY, YZ, ZZ.
 *
 * The _comps forms take ncomp separate arrays and interleave them, which is
 * the usual shape of solver data; the plain forms take one buffer already
 * laid out as (n, ncomp).
 */
h5c_status_t h5c_viz_write_point_data(h5c_viz_t *viz, const char *name,
                                      const void *buf, h5c_type_t type,
                                      size_t ncomp);

h5c_status_t h5c_viz_write_point_data_comps(h5c_viz_t *viz, const char *name,
                                            const void *const *comps,
                                            h5c_type_t type, size_t ncomp);

h5c_status_t h5c_viz_write_cell_data(h5c_viz_t *viz, const char *name,
                                     const void *buf, h5c_type_t type,
                                     size_t ncomp);

h5c_status_t h5c_viz_write_cell_data_comps(h5c_viz_t *viz, const char *name,
                                           const void *const *comps,
                                           h5c_type_t type, size_t ncomp);

/*
 * XDMF's name for a component count: "Scalar", "Vector", "Tensor6",
 * "Tensor", or NULL when XDMF has no name for it. Exposed because a caller
 * writing its own metadata needs the same mapping.
 */
const char *h5c_viz_attribute_type(size_t ncomp);

#ifdef __cplusplus
}
#endif

#endif /* H5C_VIZ_H */
