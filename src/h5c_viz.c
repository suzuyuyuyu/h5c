/*
 * Visualization HDF5 writer (scheme_version 1). Compiled only with
 * H5C_ENABLE_PARALLEL, like h5c_parallel.c.
 *
 * The file layout, what each rank contributes and the collective discipline
 * are documented in include/h5c/h5c_viz.h. h5fortran's
 * docs/USAGE-visualization.md is the format specification and
 * src/parallel/h5fort_parallel_visualization.F90 the reference implementation.
 *
 * WHY THE HELPERS BELOW LOOK LIKE h5c_parallel.c
 * ----------------------------------------------
 * agree(), make_dxpl(), the empty-selection rule and the tile plan are the
 * same devices h5c_parallel.c uses, and they solve the same problems here.
 * They are re-expressed rather than called because they are static to that
 * translation unit AND because the two layers store different things:
 * h5c_pwrite() writes a GROUP holding "data" plus "__partition__", while this
 * writer must produce PLAIN datasets at fixed paths so that h5xdmf and
 * h5fortran read the very same file. What genuinely can be shared is shared:
 * the interleave arithmetic comes from h5c__tile_rows() / h5c__pack_tile()
 * and the attribute code from h5c_attribute.c, through a borrowed
 * h5c_file_t wrapper over this file's id.
 */
#include "h5c_internal.h"

#include "h5c/h5c_viz.h"

#include <stdlib.h>
#include <string.h>

/* Fixed names of the subgroups and the geometry datasets. */
#define GEOM_GROUP  "geometry"
#define PDATA_GROUP "point_data"
#define CDATA_GROUP "cell_data"
#define NODES_NAME  "nodes"
#define CONN_NAME   "connectivity"

/* Number of coordinates per node. The scheme is always three-dimensional. */
#define NODE_COMPS 3

struct h5c_viz {
    hid_t        fid;
    h5c_file_t  *wrap;      /* borrowed wrapper over fid, for the attr code */
    h5c_status_t sticky;    /* first non-OK status seen on this writer */
    MPI_Comm     comm;      /* borrowed from the caller, valid until close */
    int          me;
    int          nprocs;

    /* --- the ONE current mesh; see close_mesh() ------------------- */
    int            have_mesh;
    char          *name;    /* owned */
    h5c_viz_kind_t kind;
    int            npe;     /* nodes per element; 1 for POLYDATA */
    size_t         num_points;
    size_t         num_cells;
    size_t         point_offset;
    size_t         cell_offset;
    size_t         total_points;
    size_t         total_cells;
    hid_t          gid_mesh;
    hid_t          gid_geom;
    hid_t          gid_pdata;
    hid_t          gid_cdata;  /* H5I_INVALID_HID for POLYDATA */
};

/* ------------------------------------------------------------------ */
/* collective agreement                                                */
/* ------------------------------------------------------------------ */

/*
 * Folds every rank's status into one, so that a rank which fails validation
 * cannot return early while the others walk into a collective HDF5 call. The
 * enum is append-only with H5C_OK == 0, so MAX picks a real failure over
 * success. Identical in spirit to agree() in h5c_parallel.c.
 */
static h5c_status_t agree(MPI_Comm comm, h5c_status_t local)
{
    int mine, worst;

    mine  = (int)local;
    worst = mine;
    if (MPI_Allreduce(&mine, &worst, 1, MPI_INT, MPI_MAX, comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI, "MPI_Allreduce failed agreeing on status");
    }
    if (worst == (int)H5C_OK) {
        return H5C_OK;
    }
    if (local != H5C_OK) {
        return local; /* keep this rank's own, more specific message */
    }
    return h5c__fail((h5c_status_t)worst,
                     "another rank reported '%s'; failing collectively",
                     h5c_status_string((h5c_status_t)worst));
}

/* Folds `status` into the writer's sticky error and returns it unchanged. */
static h5c_status_t record(h5c_viz_t *viz, h5c_status_t status)
{
    if (viz != NULL && viz->sticky == H5C_OK && status != H5C_OK) {
        viz->sticky = status;
    }
    return status;
}

/* Validates the handle. Nothing can be agreed without a communicator. */
static h5c_status_t viz_check(const h5c_viz_t *viz)
{
    if (viz == NULL || viz->fid < 0) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "visualization handle is NULL or closed");
    }
    return H5C_OK;
}

/* Validates that a mesh is current. */
static h5c_status_t mesh_check(const h5c_viz_t *viz)
{
    if (!viz->have_mesh) {
        return h5c__fail(H5C_ERR_STATE,
                         "no current mesh; call h5c_viz_begin_mesh() first");
    }
    return H5C_OK;
}

/* ------------------------------------------------------------------ */
/* paths                                                               */
/* ------------------------------------------------------------------ */

/*
 * Builds "/<mesh>[/<sub>[/<leaf>]]" for the attribute calls, which address
 * objects by path. Returns NULL on allocation failure; the caller frees.
 */
static char *mesh_path(const h5c_viz_t *viz, const char *sub, const char *leaf)
{
    size_t n;
    char  *p;

    n = 1 + strlen(viz->name) + 1;
    if (sub != NULL) {
        n += strlen(sub) + 1;
    }
    if (leaf != NULL) {
        n += strlen(leaf) + 1;
    }
    p = (char *)malloc(n);
    if (p == NULL) {
        h5c__fail(H5C_ERR_NOMEM, "cannot build the object path for '%s'",
                  viz->name);
        return NULL;
    }
    if (leaf != NULL) {
        snprintf(p, n, "/%s/%s/%s", viz->name, sub, leaf);
    } else if (sub != NULL) {
        snprintf(p, n, "/%s/%s", viz->name, sub);
    } else {
        snprintf(p, n, "/%s", viz->name);
    }
    return p;
}

/* ------------------------------------------------------------------ */
/* dataspaces and transfers                                            */
/* ------------------------------------------------------------------ */

/* Transfer property list. Visualization output is always collective. */
static hid_t make_dxpl(void)
{
    hid_t xfer;

    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (xfer < 0) {
        h5c__fail_hdf5((long)xfer, "H5Pcreate(H5P_DATASET_XFER) failed");
        return H5I_INVALID_HID;
    }
    if (H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE) < 0) {
        H5Pclose(xfer);
        h5c__fail_hdf5(-1, "H5Pset_dxpl_mpio failed");
        return H5I_INVALID_HID;
    }
    return xfer;
}

/*
 * Selects file rows [offset, offset + rows) of a (total[, ncols]) dataset and
 * builds a matching contiguous memory space.
 *
 * `rows == 0` selects nothing at all on BOTH spaces rather than a zero-length
 * hyperslab, exactly as select_block() in h5c_parallel.c does: a rank that
 * owns nothing still enters every collective call.
 */
static h5c_status_t select_rows(hid_t fsid, int drank, size_t offset,
                                size_t rows, size_t ncols, hid_t *msid_out)
{
    hsize_t start[2], count[2], mdims[2];
    hid_t   msid;

    *msid_out = H5I_INVALID_HID;

    mdims[0] = (rows > 0) ? (hsize_t)rows : 1;
    mdims[1] = (hsize_t)ncols;
    msid = H5Screate_simple(drank, mdims, NULL);
    if (msid < 0) {
        return h5c__fail_hdf5((long)msid, "cannot build the memory dataspace");
    }

    if (rows == 0) {
        if (H5Sselect_none(fsid) < 0 || H5Sselect_none(msid) < 0) {
            H5Sclose(msid);
            return h5c__fail_hdf5(-1, "H5Sselect_none failed");
        }
    } else {
        start[0] = (hsize_t)offset;
        start[1] = 0;
        count[0] = (hsize_t)rows;
        count[1] = (hsize_t)ncols;
        if (H5Sselect_hyperslab(fsid, H5S_SELECT_SET, start, NULL,
                                count, NULL) < 0) {
            H5Sclose(msid);
            return h5c__fail_hdf5(-1, "H5Sselect_hyperslab failed");
        }
    }
    *msid_out = msid;
    return H5C_OK;
}

/*
 * Creates <gid>/<name> with the file shape (total[, ncols]). On success both
 * ids are open and owned by the caller.
 */
static h5c_status_t create_dataset(hid_t gid, const char *name,
                                   h5c_type_t type, size_t total,
                                   int drank, size_t ncols,
                                   hid_t *did_out, hid_t *fsid_out)
{
    hsize_t fdims[2];
    hid_t   ftype, fsid, did;
    htri_t  present;

    *did_out  = H5I_INVALID_HID;
    *fsid_out = H5I_INVALID_HID;

    ftype = h5c__file_type(type);
    if (ftype == H5I_INVALID_HID) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "type %d has no numeric mapping for '%s'",
                         (int)type, name);
    }
    present = H5Lexists(gid, name, H5P_DEFAULT);
    if (present < 0) {
        return h5c__fail_hdf5((long)present, "H5Lexists failed for '%s'", name);
    }
    if (present > 0) {
        return h5c__fail(H5C_ERR_EXISTS, "'%s' has already been written", name);
    }

    fdims[0] = (hsize_t)total;
    fdims[1] = (hsize_t)ncols;
    fsid = H5Screate_simple(drank, fdims, NULL);
    if (fsid < 0) {
        return h5c__fail_hdf5((long)fsid,
                              "cannot build the file dataspace for '%s'", name);
    }
    did = H5Dcreate2(gid, name, ftype, fsid,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (did < 0) {
        H5Sclose(fsid);
        return h5c__fail_hdf5((long)did, "cannot create the dataset '%s'",
                              name);
    }
    *did_out  = did;
    *fsid_out = fsid;
    return H5C_OK;
}

/*
 * Writes this rank's rows of an already-contiguous buffer in one collective
 * transfer. `buf` may be NULL when `rows == 0`.
 */
static h5c_status_t write_block(hid_t gid, const char *name, const void *buf,
                                h5c_type_t type, size_t rows, size_t offset,
                                size_t total, int drank, size_t ncols)
{
    h5c_status_t st;
    hid_t        did, fsid, msid = H5I_INVALID_HID, xfer;
    char         dummy = 0;

    if ((st = create_dataset(gid, name, type, total, drank, ncols,
                             &did, &fsid)) != H5C_OK) {
        return st;
    }
    xfer = make_dxpl();
    if (xfer == H5I_INVALID_HID) {
        H5Dclose(did);
        H5Sclose(fsid);
        return H5C_ERR_HDF5;
    }
    if ((st = select_rows(fsid, drank, offset, rows, ncols, &msid)) == H5C_OK) {
        if (H5Dwrite(did, h5c__mem_type(type), msid, fsid, xfer,
                     (buf != NULL) ? buf : (const void *)&dummy) < 0) {
            st = h5c__fail_hdf5(-1, "H5Dwrite failed for '%s'", name);
        }
        H5Sclose(msid);
    }
    H5Pclose(xfer);
    H5Dclose(did);
    H5Sclose(fsid);
    return st;
}

/* ------------------------------------------------------------------ */
/* staged (row-tiled) writes                                           */
/* ------------------------------------------------------------------ */

/*
 * Rows that are not contiguous in the caller's memory are staged tile by
 * tile, so the staging buffer stays bounded by h5c_pack_limit() while every
 * transfer remains contiguous in the file. `stage` gathers rows
 * [row0, row0 + rows) of the caller's data into `dst`.
 */
typedef void (*stage_fn)(char *dst, size_t row0, size_t rows, const void *ctx);

/* Component interleaving; ctx is a comps_ctx. Delegates to h5c__pack_tile. */
typedef struct {
    const void *const *comps;
    size_t             ncomp;
    size_t             esize;
} comps_ctx;

static void stage_comps(char *dst, size_t row0, size_t rows, const void *ctx)
{
    const comps_ctx *c = (const comps_ctx *)ctx;

    h5c__pack_tile(dst, c->comps, c->ncomp, row0, rows, c->esize);
}

/*
 * Connectivity: copy plus this rank's node offset, in the caller's own
 * element type. The offset addition happens HERE and never in the caller's
 * buffer, which stays const throughout.
 */
typedef struct {
    const void *conn;
    size_t      npe;
    size_t      offset;
    h5c_type_t  type;
} conn_ctx;

/* Reads element `i` of a connectivity buffer of `type` as int64. */
static int64_t conn_get(const void *conn, size_t i, h5c_type_t type)
{
    switch (type) {
    case H5C_I8:  return (int64_t)((const int8_t  *)conn)[i];
    case H5C_I16: return (int64_t)((const int16_t *)conn)[i];
    case H5C_I32: return (int64_t)((const int32_t *)conn)[i];
    default:      return          ((const int64_t *)conn)[i];
    }
}

/* Stores `v` as element `i` of a connectivity buffer of `type`. */
static void conn_put(void *dst, size_t i, h5c_type_t type, int64_t v)
{
    switch (type) {
    case H5C_I8:  ((int8_t  *)dst)[i] = (int8_t)v;  break;
    case H5C_I16: ((int16_t *)dst)[i] = (int16_t)v; break;
    case H5C_I32: ((int32_t *)dst)[i] = (int32_t)v; break;
    default:      ((int64_t *)dst)[i] = v;          break;
    }
}

/* Largest value the connectivity element type can hold. */
static int64_t conn_max(h5c_type_t type)
{
    switch (type) {
    case H5C_I8:  return 127;
    case H5C_I16: return 32767;
    case H5C_I32: return 2147483647;
    default:      return 9223372036854775807LL;
    }
}

static void stage_conn(char *dst, size_t row0, size_t rows, const void *ctx)
{
    const conn_ctx *c = (const conn_ctx *)ctx;
    size_t          first = row0 * c->npe;
    size_t          n     = rows * c->npe;
    size_t          i;

    for (i = 0; i < n; i++) {
        conn_put(dst, i, c->type,
                 conn_get(c->conn, first + i, c->type) + (int64_t)c->offset);
    }
}

/*
 * The tile count must be IDENTICAL on every rank: local row counts differ, so
 * a locally derived count would make ranks issue different numbers of
 * collective transfers and deadlock. Ranks that run out of rows early still
 * enter every call, with an empty selection. Same device as agree_tiles() in
 * h5c_parallel.c; the floor of 1 keeps one transfer when no rank owns rows.
 */
static h5c_status_t agree_tiles(MPI_Comm comm, size_t n, size_t rows,
                                long long *ntiles)
{
    long long mine, most;

    mine = (rows > 0) ? (long long)((n + rows - 1) / rows) : 0;
    most = mine;
    if (MPI_Allreduce(&mine, &most, 1, MPI_LONG_LONG, MPI_MAX,
                      comm) != MPI_SUCCESS) {
        *ntiles = 1;
        return h5c__fail(H5C_ERR_MPI,
                         "MPI_Allreduce failed agreeing on the tile count");
    }
    *ntiles = (most < 1) ? 1 : most;
    return H5C_OK;
}

/*
 * Creates <gid>/<name> and writes this rank's rows through `stage`, one
 * collective transfer per tile. A transfer error is remembered but does NOT
 * leave the loop: every rank must issue the same number of collective calls.
 *
 * The returned status is agreed across the communicator, so all ranks leave
 * with the same verdict.
 */
static h5c_status_t write_staged(MPI_Comm comm, hid_t gid, const char *name,
                                 h5c_type_t type, size_t rows_local,
                                 size_t offset, size_t total, size_t ncols,
                                 stage_fn stage, const void *ctx)
{
    h5c_status_t st = H5C_OK;
    hid_t        did = H5I_INVALID_HID, fsid = H5I_INVALID_HID;
    hid_t        xfer = H5I_INVALID_HID;
    char        *buf = NULL;
    size_t       esize, row_bytes, rows;
    long long    ntiles = 1, t;
    int          drank;

    esize     = h5c_type_size(type);
    drank     = (ncols > 1) ? 2 : 1;
    row_bytes = (esize > 0) ? ncols * esize : 1;
    rows      = h5c__tile_rows(rows_local, row_bytes);
    if (esize == 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "type %d cannot be staged for '%s'", (int)type, name);
    } else if (rows > 0) {
        buf = (char *)malloc(rows * row_bytes);
        if (buf == NULL) {
            st = h5c__fail(H5C_ERR_NOMEM,
                           "cannot allocate %lu bytes to stage '%s'",
                           (unsigned long)(rows * row_bytes), name);
        }
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        free(buf);
        return st;
    }
    if ((st = agree_tiles(comm, rows_local, rows, &ntiles)) != H5C_OK) {
        free(buf);
        return st;
    }

    st = create_dataset(gid, name, type, total, drank, ncols, &did, &fsid);
    if (st == H5C_OK) {
        xfer = make_dxpl();
        if (xfer == H5I_INVALID_HID) {
            st = H5C_ERR_HDF5;
        }
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        if (xfer >= 0) { H5Pclose(xfer); }
        if (did  >= 0) { H5Dclose(did);  }
        if (fsid >= 0) { H5Sclose(fsid); }
        free(buf);
        return st;
    }

    for (t = 0; t < ntiles; t++) {
        size_t row0 = (size_t)t * rows;
        size_t take = 0;
        char   dummy = 0;
        hid_t  msid;

        if (rows > 0 && row0 < rows_local) {
            take = (rows_local - row0 < rows) ? rows_local - row0 : rows;
        }
        if (select_rows(fsid, drank, offset + row0, take, ncols,
                        &msid) != H5C_OK) {
            /* Unreachable in practice; nothing is left to select. */
            if (st == H5C_OK) {
                st = H5C_ERR_HDF5;
            }
            break;
        }
        if (take > 0) {
            stage(buf, row0, take, ctx);
        }
        if (H5Dwrite(did, h5c__mem_type(type), msid, fsid, xfer,
                     (buf != NULL) ? (const void *)buf
                                   : (const void *)&dummy) < 0 &&
            st == H5C_OK) {
            st = h5c__fail_hdf5(-1, "H5Dwrite failed for tile %lld of '%s'",
                                t, name);
        }
        H5Sclose(msid);
    }

    if (xfer >= 0) { H5Pclose(xfer); }
    if (did  >= 0) { H5Dclose(did);  }
    if (fsid >= 0) { H5Sclose(fsid); }
    free(buf);
    return agree(comm, st);
}

/* ------------------------------------------------------------------ */
/* attribute_type                                                      */
/* ------------------------------------------------------------------ */

const char *h5c_viz_attribute_type(size_t ncomp)
{
    switch (ncomp) {
    case 1: return "Scalar";
    case 3: return "Vector";
    case 6: return "Tensor6";
    case 9: return "Tensor";
    default: return NULL;
    }
}

/* ------------------------------------------------------------------ */
/* lifecycle                                                           */
/* ------------------------------------------------------------------ */

h5c_status_t h5c_viz_open(const char *path, double time,
                          MPI_Comm comm, MPI_Info info, h5c_viz_t **out)
{
    h5c_status_t st;
    h5c_viz_t   *viz;
    hid_t        fapl, fid;
    int32_t      scheme = (int32_t)H5C_SCHEME_VERSION;

    if (out == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_viz_open: out is NULL");
    }
    *out = NULL;
    if (path == NULL || path[0] == '\0') {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_viz_open: empty path");
    }
    if (comm == MPI_COMM_NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_viz_open: MPI_COMM_NULL");
    }
    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }

    fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl < 0) {
        return h5c__fail_hdf5((long)fapl, "H5Pcreate(H5P_FILE_ACCESS) failed");
    }
    /*
     * DELIBERATE, and copied from h5fortran: no collective metadata
     * properties. Every rank writes every attribute symmetrically instead,
     * which is what keeps the metadata consistent; mixing rank-0-only
     * attribute writes with collective metadata operations corrupts the
     * metadata checksums.
     */
    if (H5Pset_fapl_mpio(fapl, comm, info) < 0) {
        H5Pclose(fapl);
        return h5c__fail_hdf5(-1, "H5Pset_fapl_mpio failed for '%s'", path);
    }
    fid = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    if (fid < 0) {
        return h5c__fail_hdf5((long)fid, "cannot create '%s' in parallel",
                              path);
    }

    viz = (h5c_viz_t *)calloc(1, sizeof *viz);
    if (viz == NULL) {
        H5Fclose(fid);
        return h5c__fail(H5C_ERR_NOMEM, "h5c_viz_open: allocation failed");
    }
    viz->fid       = fid;
    viz->sticky    = H5C_OK;
    viz->comm      = comm;   /* borrowed; valid until h5c_viz_close() */
    viz->gid_mesh  = H5I_INVALID_HID;
    viz->gid_geom  = H5I_INVALID_HID;
    viz->gid_pdata = H5I_INVALID_HID;
    viz->gid_cdata = H5I_INVALID_HID;

    /* The wrapper is borrowed: closing it never closes `fid`. */
    if ((st = h5c_file_from_hid(fid, &viz->wrap)) != H5C_OK) {
        H5Fclose(fid);
        free(viz);
        return st;
    }
    if (MPI_Comm_rank(comm, &viz->me) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &viz->nprocs) != MPI_SUCCESS) {
        h5c_close(viz->wrap);
        H5Fclose(fid);
        free(viz);
        return h5c__fail(H5C_ERR_MPI, "MPI_Comm_rank/size failed");
    }

    /* Root metadata, written by every rank so that all ranks agree. */
    st = h5c_write_attr_scalar(viz->wrap, "/", "scheme_version",
                               &scheme, H5C_I32);
    if (st == H5C_OK) {
        st = h5c_write_attr_scalar(viz->wrap, "/", "time", &time, H5C_F64);
    }
    if ((st = agree(comm, st)) != H5C_OK) {
        h5c_close(viz->wrap);
        H5Fclose(fid);
        free(viz);
        return st;
    }

    *out = viz;
    return H5C_OK;
}

/* Releases the current mesh's group ids. Leaves no dangling id behind. */
static void close_mesh(h5c_viz_t *viz)
{
    if (viz->gid_cdata >= 0) { H5Gclose(viz->gid_cdata); }
    if (viz->gid_pdata >= 0) { H5Gclose(viz->gid_pdata); }
    if (viz->gid_geom  >= 0) { H5Gclose(viz->gid_geom);  }
    if (viz->gid_mesh  >= 0) { H5Gclose(viz->gid_mesh);  }
    viz->gid_cdata = H5I_INVALID_HID;
    viz->gid_pdata = H5I_INVALID_HID;
    viz->gid_geom  = H5I_INVALID_HID;
    viz->gid_mesh  = H5I_INVALID_HID;
    free(viz->name);
    viz->name      = NULL;
    viz->have_mesh = 0;
}

h5c_status_t h5c_viz_close(h5c_viz_t *viz)
{
    h5c_status_t st = H5C_OK;

    if (viz == NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_viz_close: viz is NULL");
    }
    close_mesh(viz);
    if (viz->fid >= 0) {
        if (H5Fflush(viz->fid, H5F_SCOPE_GLOBAL) < 0) {
            st = h5c__fail_hdf5(-1, "H5Fflush failed");
        }
        if (H5Fclose(viz->fid) < 0 && st == H5C_OK) {
            st = h5c__fail_hdf5(-1, "H5Fclose failed");
        }
    }
    if (viz->wrap != NULL) {
        h5c_close(viz->wrap);  /* borrowed wrapper; frees only the wrapper */
    }
    free(viz);
    return st;
}

h5c_status_t h5c_viz_status(const h5c_viz_t *viz)
{
    return (viz == NULL) ? H5C_ERR_INVALID_ARG : viz->sticky;
}

/* ------------------------------------------------------------------ */
/* meshes                                                              */
/* ------------------------------------------------------------------ */

/* Opens <loc>/<name>, creating it when absent. */
static h5c_status_t open_or_create_group(hid_t loc, const char *name,
                                         hid_t *out)
{
    htri_t present;
    hid_t  gid;

    *out    = H5I_INVALID_HID;
    present = H5Lexists(loc, name, H5P_DEFAULT);
    if (present < 0) {
        return h5c__fail_hdf5((long)present, "H5Lexists failed for '%s'", name);
    }
    if (present > 0) {
        gid = H5Gopen2(loc, name, H5P_DEFAULT);
        if (gid < 0) {
            return h5c__fail_hdf5((long)gid, "cannot open the group '%s'",
                                  name);
        }
    } else {
        gid = H5Gcreate2(loc, name, h5c__lcpl(), H5P_DEFAULT, H5P_DEFAULT);
        if (gid < 0) {
            return h5c__fail_hdf5((long)gid, "cannot create the group '%s'",
                                  name);
        }
    }
    *out = gid;
    return H5C_OK;
}

/*
 * Derives this rank's offsets and the totals from the local counts.
 *
 * DELIBERATE: h5fortran MPI_Allgathers both counts and then sums the ranks
 * below its own. The prefix sum over the ranks in order is exactly what
 * MPI_Exscan computes, and MPI_Allreduce the totals, so the same numbers come
 * out of two allocation-free calls. That matters here beyond tidiness: a
 * calloc for nprocs entries could fail on one rank only, and a rank that
 * skipped the gather while the others entered it would hang.
 *
 * Point and cell offsets are independent: both are carried in one pair.
 */
static h5c_status_t gather_counts(h5c_viz_t *viz, size_t np, size_t nc)
{
    int64_t mine[2], before[2], total[2];

    mine[0]   = (int64_t)np;
    mine[1]   = (int64_t)nc;
    before[0] = 0;
    before[1] = 0;  /* MPI_Exscan leaves rank 0's result undefined */

    if (MPI_Exscan(mine, before, 2, MPI_INT64_T, MPI_SUM,
                   viz->comm) != MPI_SUCCESS ||
        MPI_Allreduce(mine, total, 2, MPI_INT64_T, MPI_SUM,
                      viz->comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI,
                         "MPI_Exscan/Allreduce failed collecting the counts");
    }
    if (viz->me == 0) {
        before[0] = 0;
        before[1] = 0;
    }

    viz->point_offset = (size_t)before[0];
    viz->cell_offset  = (size_t)before[1];
    viz->total_points = (size_t)total[0];
    viz->total_cells  = (size_t)total[1];
    return H5C_OK;
}

h5c_status_t h5c_viz_begin_mesh(h5c_viz_t *viz, const h5c_viz_mesh_t *mesh)
{
    h5c_status_t st;
    const char  *name, *topology;
    char        *path = NULL;
    int          npe, created;
    int32_t      npe_attr;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return record(viz, st);
    }
    if ((st = viz_check(viz)) != H5C_OK) {
        return st;  /* no communicator, so nothing can be agreed */
    }

    /*
     * The previous mesh is deselected FIRST, before anything can fail:
     * exactly one mesh is ever current, and a caller must not be left
     * writing into the previous mesh because its begin_mesh() was rejected.
     * Every failure below therefore leaves no current mesh, and no group id
     * dangling either.
     */
    close_mesh(viz);

    /* --- local validation, agreed before any HDF5 or MPI call ------ */
    st = H5C_OK;
    name     = H5C_VIZ_DEFAULT_UGRID_NAME;
    topology = H5C_VIZ_DEFAULT_TOPOLOGY;
    npe      = H5C_VIZ_DEFAULT_NODES_PER_ELEM;
    if (mesh == NULL) {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "h5c_viz_begin_mesh: mesh is NULL");
    } else if (mesh->kind == H5C_VIZ_UNSTRUCTURED) {
        if (mesh->name != NULL && mesh->name[0] != '\0') {
            name = mesh->name;
        }
        if (mesh->topology != NULL && mesh->topology[0] != '\0') {
            topology = mesh->topology;
        }
        if (mesh->nodes_per_element > 0) {
            npe = mesh->nodes_per_element;
        } else if (mesh->nodes_per_element < 0) {
            st = h5c__fail(H5C_ERR_INVALID_ARG,
                           "nodes_per_element %d is negative",
                           mesh->nodes_per_element);
        }
    } else if (mesh->kind == H5C_VIZ_POLYDATA) {
        /*
         * A point cloud has no elements: h5fortran forces Polyvertex and 1
         * whatever the caller asked for, and the same values must appear here
         * or h5xdmf would emit a topology the file cannot back.
         */
        name     = H5C_VIZ_DEFAULT_POLYDATA_NAME;
        topology = "Polyvertex";
        npe      = 1;
        if (mesh->name != NULL && mesh->name[0] != '\0') {
            name = mesh->name;
        }
        if (mesh->num_cells != 0) {
            st = h5c__fail(H5C_ERR_INVALID_ARG,
                           "POLYDATA mesh '%s' has num_cells %lu; "
                           "a point cloud owns no cells",
                           name, (unsigned long)mesh->num_cells);
        }
    } else {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "unknown mesh kind %d",
                       (int)(mesh == NULL ? 0 : mesh->kind));
    }
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        return record(viz, st);
    }

    viz->kind       = mesh->kind;
    viz->npe        = npe;
    viz->num_points = mesh->num_points;
    viz->num_cells  = mesh->num_cells;
    viz->name       = (char *)malloc(strlen(name) + 1);
    if (viz->name == NULL) {
        st = h5c__fail(H5C_ERR_NOMEM, "cannot copy the mesh name '%s'", name);
    } else {
        memcpy(viz->name, name, strlen(name) + 1);
    }
    /* Agreed before gather_counts(), which every rank must enter together. */
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        close_mesh(viz);
        return record(viz, st);
    }
    if ((st = gather_counts(viz, mesh->num_points, mesh->num_cells)) != H5C_OK) {
        close_mesh(viz);
        return record(viz, st);
    }

    /* --- groups, collective ---------------------------------------- */
    {
        htri_t present = H5Lexists(viz->fid, viz->name, H5P_DEFAULT);

        if (present < 0) {
            st = h5c__fail_hdf5((long)present, "H5Lexists failed for '%s'",
                                viz->name);
        }
        created = (present == 0);
    }
    if (st == H5C_OK) {
        st = open_or_create_group(viz->fid, viz->name, &viz->gid_mesh);
    }
    if (st == H5C_OK) {
        st = open_or_create_group(viz->gid_mesh, GEOM_GROUP, &viz->gid_geom);
    }
    if (st == H5C_OK) {
        st = open_or_create_group(viz->gid_mesh, PDATA_GROUP, &viz->gid_pdata);
    }
    if (st == H5C_OK && viz->kind == H5C_VIZ_UNSTRUCTURED) {
        st = open_or_create_group(viz->gid_mesh, CDATA_GROUP, &viz->gid_cdata);
    }

    /*
     * Topology metadata is stamped only when this call created the group, so
     * that re-selecting a mesh does not rewrite (delete and recreate) the
     * attributes under every rank.
     */
    if (st == H5C_OK && created) {
        path = mesh_path(viz, NULL, NULL);
        if (path == NULL) {
            st = H5C_ERR_NOMEM;
        } else {
            npe_attr = (int32_t)npe;
            st = h5c_write_attr_str(viz->wrap, path, "topology_type", topology);
            if (st == H5C_OK) {
                st = h5c_write_attr_scalar(viz->wrap, path,
                                           "nodes_per_element", &npe_attr,
                                           H5C_I32);
            }
            free(path);
        }
    }
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        close_mesh(viz);
        return record(viz, st);
    }

    viz->have_mesh = 1;
    return H5C_OK;
}

h5c_status_t h5c_viz_offsets(const h5c_viz_t *viz,
                             size_t *point_offset, size_t *cell_offset)
{
    h5c_status_t st;

    if ((st = viz_check(viz)) != H5C_OK) {
        return st;
    }
    if ((st = mesh_check(viz)) != H5C_OK) {
        return st;
    }
    if (point_offset != NULL) {
        *point_offset = viz->point_offset;
    }
    if (cell_offset != NULL) {
        *cell_offset = viz->cell_offset;
    }
    return H5C_OK;
}

/* ------------------------------------------------------------------ */
/* geometry                                                            */
/* ------------------------------------------------------------------ */

/* Shared entry checks for every write call below. */
static h5c_status_t begin_write(h5c_viz_t *viz)
{
    h5c_status_t st;

    if ((st = h5c__ensure_init()) != H5C_OK) {
        return st;
    }
    if ((st = viz_check(viz)) != H5C_OK) {
        return st;
    }
    return mesh_check(viz);
}

/* F32 or F64: coordinates are real by definition. */
static h5c_status_t check_coord_type(h5c_type_t type)
{
    if (type != H5C_F32 && type != H5C_F64) {
        return h5c__fail(H5C_ERR_INVALID_ARG,
                         "nodes must be H5C_F32 or H5C_F64, not type %d",
                         (int)type);
    }
    return H5C_OK;
}

h5c_status_t h5c_viz_write_nodes(h5c_viz_t *viz, const void *nodes,
                                 h5c_type_t type)
{
    h5c_status_t st;

    if ((st = begin_write(viz)) != H5C_OK) {
        return record(viz, st);
    }
    if ((st = check_coord_type(type)) == H5C_OK &&
        nodes == NULL && viz->num_points > 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "nodes is NULL but this rank owns %lu points",
                       (unsigned long)viz->num_points);
    }
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        return record(viz, st);
    }

    st = write_block(viz->gid_geom, NODES_NAME, nodes, type, viz->num_points,
                     viz->point_offset, viz->total_points, 2, NODE_COMPS);
    return record(viz, agree(viz->comm, st));
}

h5c_status_t h5c_viz_write_nodes_comps(h5c_viz_t *viz,
                                       const void *const *xyz,
                                       h5c_type_t type)
{
    h5c_status_t st;
    comps_ctx    ctx;
    size_t       c;

    if ((st = begin_write(viz)) != H5C_OK) {
        return record(viz, st);
    }
    st = check_coord_type(type);
    if (st == H5C_OK && xyz == NULL) {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "xyz is NULL");
    }
    if (st == H5C_OK && viz->num_points > 0) {
        for (c = 0; c < NODE_COMPS; c++) {
            if (xyz[c] == NULL) {
                st = h5c__fail(H5C_ERR_INVALID_ARG, "xyz[%lu] is NULL",
                               (unsigned long)c);
                break;
            }
        }
    }
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        return record(viz, st);
    }

    ctx.comps = xyz;
    ctx.ncomp = NODE_COMPS;
    ctx.esize = h5c_type_size(type);
    st = write_staged(viz->comm, viz->gid_geom, NODES_NAME, type,
                      viz->num_points, viz->point_offset, viz->total_points,
                      NODE_COMPS, stage_comps, &ctx);
    return record(viz, st);
}

h5c_status_t h5c_viz_write_connectivity(h5c_viz_t *viz, const void *conn,
                                        h5c_type_t type)
{
    h5c_status_t st;
    conn_ctx     ctx;
    size_t       n, i;
    int64_t      limit;

    if ((st = begin_write(viz)) != H5C_OK) {
        return record(viz, st);
    }

    st = H5C_OK;
    if (type != H5C_I8 && type != H5C_I16 &&
        type != H5C_I32 && type != H5C_I64) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "connectivity must be an integer type, not type %d",
                       (int)type);
    } else if (viz->kind != H5C_VIZ_UNSTRUCTURED) {
        st = h5c__fail(H5C_ERR_STATE,
                       "mesh '%s' is a point cloud and has no connectivity",
                       viz->name);
    } else if (conn == NULL && viz->num_cells > 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "conn is NULL but this rank owns %lu cells",
                       (unsigned long)viz->num_cells);
    }

    /*
     * The caller passes RANK-LOCAL 0-origin indices and the file stores
     * file-global ids, so this rank's point offset is added on the way out.
     * Validate first, in full: an index outside [0, num_points) means the
     * caller handed over 1-origin or already-global ids, and adding the
     * offset to those would silently produce a mesh pointing at another
     * rank's nodes. The check also rejects a global id too large for the
     * element type the caller chose, which is the one way the addition
     * itself can go wrong (h5fortran wraps around instead).
     */
    if (st == H5C_OK && viz->num_cells > 0) {
        n     = viz->num_cells * (size_t)viz->npe;
        limit = conn_max(type);
        for (i = 0; i < n; i++) {
            int64_t v = conn_get(conn, i, type);

            if (v < 0 || (uint64_t)v >= (uint64_t)viz->num_points) {
                st = h5c__fail(H5C_ERR_INVALID_ARG,
                               "connectivity[%lu] is %lld, outside "
                               "[0, %lu): rank-local 0-origin node indices "
                               "are required",
                               (unsigned long)i, (long long)v,
                               (unsigned long)viz->num_points);
                break;
            }
            if (v > limit - (int64_t)viz->point_offset) {
                st = h5c__fail(H5C_ERR_INVALID_ARG,
                               "global node id %lld does not fit the "
                               "connectivity element type (max %lld); "
                               "use a wider integer type",
                               (long long)(v + (int64_t)viz->point_offset),
                               (long long)limit);
                break;
            }
        }
    }
    /* Agreed BEFORE the first collective HDF5 call, so nobody deadlocks. */
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        return record(viz, st);
    }

    ctx.conn   = conn;
    ctx.npe    = (size_t)viz->npe;
    ctx.offset = viz->point_offset;
    ctx.type   = type;
    st = write_staged(viz->comm, viz->gid_geom, CONN_NAME, type,
                      viz->num_cells, viz->cell_offset, viz->total_cells,
                      (size_t)viz->npe, stage_conn, &ctx);
    return record(viz, st);
}

/* ------------------------------------------------------------------ */
/* fields                                                              */
/* ------------------------------------------------------------------ */

/* Stamps attribute_type on a field, for the component counts XDMF names. */
static h5c_status_t write_field_attr(h5c_viz_t *viz, const char *group,
                                     const char *name, size_t ncomp)
{
    const char  *kind = h5c_viz_attribute_type(ncomp);
    h5c_status_t st;
    char        *path;

    if (kind == NULL) {
        return H5C_OK;  /* XDMF has no name for this count; see h5c_viz.h */
    }
    path = mesh_path(viz, group, name);
    if (path == NULL) {
        return H5C_ERR_NOMEM;
    }
    st = h5c_write_attr_str(viz->wrap, path, "attribute_type", kind);
    free(path);
    return st;
}

/*
 * Common half of the four field entry points. `use_comps` selects the staged
 * path (and then `comps` is the caller's array) from the contiguous one (and
 * then `buf` is the caller's buffer).
 *
 * The flag is what decides, NOT `comps != NULL`: a NULL array must be
 * rejected on viz->comm like any other bad argument, and a rank that guessed
 * the path from the pointer would take a different branch from its peers.
 */
static h5c_status_t write_field(h5c_viz_t *viz, int cell_data,
                                const char *name, const void *buf,
                                const void *const *comps, int use_comps,
                                h5c_type_t type, size_t ncomp)
{
    h5c_status_t st;
    comps_ctx    ctx;
    const char  *group;
    hid_t        gid;
    size_t       rows, offset, total, c;

    if ((st = begin_write(viz)) != H5C_OK) {
        return record(viz, st);
    }

    group  = cell_data ? CDATA_GROUP : PDATA_GROUP;
    gid    = cell_data ? viz->gid_cdata : viz->gid_pdata;
    rows   = cell_data ? viz->num_cells : viz->num_points;
    offset = cell_data ? viz->cell_offset : viz->point_offset;
    total  = cell_data ? viz->total_cells : viz->total_points;

    st = H5C_OK;
    if (name == NULL || name[0] == '\0') {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "empty field name");
    } else if (cell_data && viz->kind != H5C_VIZ_UNSTRUCTURED) {
        st = h5c__fail(H5C_ERR_STATE,
                       "mesh '%s' is a point cloud and has no cell data",
                       viz->name);
    } else if (ncomp == 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "field '%s' has ncomp 0", name);
    } else if (h5c__file_type(type) == H5I_INVALID_HID ||
               h5c_type_size(type) == 0) {
        st = h5c__fail(H5C_ERR_INVALID_ARG,
                       "type %d has no numeric mapping for field '%s'",
                       (int)type, name);
    } else if (use_comps && comps == NULL) {
        st = h5c__fail(H5C_ERR_INVALID_ARG, "field '%s': comps is NULL", name);
    } else if (rows > 0) {
        if (!use_comps && buf == NULL) {
            st = h5c__fail(H5C_ERR_INVALID_ARG,
                           "field '%s' buffer is NULL but this rank owns "
                           "%lu rows", name, (unsigned long)rows);
        }
        for (c = 0; use_comps && c < ncomp; c++) {
            if (comps[c] == NULL) {
                st = h5c__fail(H5C_ERR_INVALID_ARG,
                               "field '%s': comps[%lu] is NULL", name,
                               (unsigned long)c);
                break;
            }
        }
    }
    if ((st = agree(viz->comm, st)) != H5C_OK) {
        return record(viz, st);
    }

    if (use_comps) {
        ctx.comps = comps;
        ctx.ncomp = ncomp;
        ctx.esize = h5c_type_size(type);
        st = write_staged(viz->comm, gid, name, type, rows, offset, total,
                          ncomp, stage_comps, &ctx);
    } else {
        st = agree(viz->comm,
                   write_block(gid, name, buf, type, rows, offset, total,
                               (ncomp > 1) ? 2 : 1, ncomp));
    }
    if (st != H5C_OK) {
        return record(viz, st);
    }

    /* Every rank writes the attribute, as for the root metadata. */
    return record(viz, agree(viz->comm,
                             write_field_attr(viz, group, name, ncomp)));
}

h5c_status_t h5c_viz_write_point_data(h5c_viz_t *viz, const char *name,
                                      const void *buf, h5c_type_t type,
                                      size_t ncomp)
{
    return write_field(viz, 0, name, buf, NULL, 0, type, ncomp);
}

h5c_status_t h5c_viz_write_point_data_comps(h5c_viz_t *viz, const char *name,
                                            const void *const *comps,
                                            h5c_type_t type, size_t ncomp)
{
    return write_field(viz, 0, name, NULL, comps, 1, type, ncomp);
}

h5c_status_t h5c_viz_write_cell_data(h5c_viz_t *viz, const char *name,
                                     const void *buf, h5c_type_t type,
                                     size_t ncomp)
{
    return write_field(viz, 1, name, buf, NULL, 0, type, ncomp);
}

h5c_status_t h5c_viz_write_cell_data_comps(h5c_viz_t *viz, const char *name,
                                           const void *const *comps,
                                           h5c_type_t type, size_t ncomp)
{
    return write_field(viz, 1, name, NULL, comps, 1, type, ncomp);
}
