#include "h5c_viz_internal.h"
#include "h5c/h5c_viz_mpi.h"

#include <stdlib.h>

typedef struct {
    MPI_Comm comm;
    int me, nprocs;
} viz_context;

static h5c_status_t agree(const void *context, h5c_status_t local)
{
    const viz_context *ctx = context;
    int mine = (int)local, worst = mine;
    if (MPI_Allreduce(&mine, &worst, 1, MPI_INT, MPI_MAX,
                      ctx->comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI, "MPI_Allreduce failed agreeing on status");
    }
    if (local != H5C_OK || worst == H5C_OK) { return local; }
    return h5c__fail((h5c_status_t)worst,
                     "another rank reported '%s'; failing collectively",
                     h5c_status_string((h5c_status_t)worst));
}

static hid_t make_dxpl(void)
{
    hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
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

static h5c_status_t select_none(hid_t fsid, hid_t msid)
{
    if (H5Sselect_none(fsid) < 0 || H5Sselect_none(msid) < 0) {
        return h5c__fail_hdf5(-1, "H5Sselect_none failed");
    }
    return H5C_OK;
}

static h5c_status_t gather_counts(const void *context, h5c_viz_t *viz,
                                   size_t np, size_t nc)
{
    const viz_context *ctx = context;
    int64_t mine[2] = { (int64_t)np, (int64_t)nc };
    int64_t offsets[2], totals[2];
    if (MPI_Exscan(mine, offsets, 2, MPI_INT64_T, MPI_SUM,
                   ctx->comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI, "MPI_Exscan failed collecting offsets");
    }
    /* MPI_Exscan leaves rank 0's result undefined. */
    if (ctx->me == 0) { offsets[0] = offsets[1] = 0; }
    if (MPI_Allreduce(mine, totals, 2, MPI_INT64_T, MPI_SUM,
                      ctx->comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI, "MPI_Allreduce failed collecting totals");
    }
    viz->point_offset = (size_t)offsets[0];
    viz->cell_offset = (size_t)offsets[1];
    viz->total_points = (size_t)totals[0];
    viz->total_cells = (size_t)totals[1];
    return H5C_OK;
}

static h5c_status_t agree_tiles(const void *context, long long *ntiles)
{
    const viz_context *ctx = context;
    long long mine = *ntiles;
    if (MPI_Allreduce(&mine, ntiles, 1, MPI_LONG_LONG, MPI_MAX,
                      ctx->comm) != MPI_SUCCESS) {
        return h5c__fail(H5C_ERR_MPI, "MPI_Allreduce failed agreeing on tile count");
    }
    return H5C_OK;
}

static const h5c_viz_ops ops = {
    agree, make_dxpl, select_none, gather_counts, agree_tiles
};

h5c_status_t h5c_viz_popen(const char *path, double time, MPI_Comm comm,
                           MPI_Info info, h5c_viz_t **out)
{
    viz_context initial = { .comm = comm }, *ctx = NULL;
    h5c_status_t st;
    hid_t fapl = H5I_INVALID_HID;

    if (out) { *out = NULL; }
    if (comm == MPI_COMM_NULL) {
        return h5c__fail(H5C_ERR_INVALID_ARG, "h5c_viz_popen: MPI_COMM_NULL");
    }
    st = h5c__ensure_init();
    if (st == H5C_OK &&
        (MPI_Comm_rank(comm, &initial.me) != MPI_SUCCESS ||
         MPI_Comm_size(comm, &initial.nprocs) != MPI_SUCCESS)) {
        st = h5c__fail(H5C_ERR_MPI, "MPI_Comm_rank/size failed");
    }
    if (st == H5C_OK) {
        ctx = malloc(sizeof *ctx);
        if (!ctx) { st = h5c__fail(H5C_ERR_NOMEM, "cannot allocate writer context"); }
        else { *ctx = initial; }
    }
    if (st == H5C_OK) {
        fapl = H5Pcreate(H5P_FILE_ACCESS);
        /* Every rank writes attributes; do not enable collective metadata. */
        if (fapl < 0 || H5Pset_fapl_mpio(fapl, comm, info) < 0) {
            st = h5c__fail_hdf5(-1, "cannot configure parallel file access");
        }
    }
    st = agree(&initial, st);
    if (st == H5C_OK) {
        st = h5c__viz_open(path, time, fapl, &ops, ctx, out);
    } else {
        free(ctx);
    }
    if (fapl >= 0) { H5Pclose(fapl); }
    return st;
}
