/* Shared writer state and backend hooks. Not installed. */
#ifndef H5C_VIZ_INTERNAL_H
#define H5C_VIZ_INTERNAL_H

#include "h5c/h5c_viz.h"
#include "h5c_internal.h"

typedef struct {
    h5c_status_t (*agree)(const void*, h5c_status_t);
    hid_t (*make_dxpl)(void);
    h5c_status_t (*select_none)(hid_t, hid_t);
    h5c_status_t (*gather_counts)(const void*, h5c_viz_t*, size_t, size_t);
    h5c_status_t (*agree_tiles)(const void*, long long*);
} h5c_viz_ops;

struct h5c_viz {
    hid_t fid;
    h5c_file_t* wrap;    /* borrowed wrapper over fid, for the attr code */
    h5c_status_t sticky; /* first non-OK status seen on this writer */
    const h5c_viz_ops* ops;
    void* context;

    /* the ONE current mesh; see close_mesh() */
    int have_mesh;
    char* name; /* owned */
    h5c_viz_kind_t kind;
    int npe; /* nodes per element; 1 for POLYDATA */
    size_t num_points;
    size_t num_cells;
    size_t point_offset;
    size_t cell_offset;
    size_t total_points;
    size_t total_cells;
    hid_t gid_mesh;
    hid_t gid_geom;
    hid_t gid_pdata;
    hid_t gid_cdata; /* H5I_INVALID_HID for POLYDATA */
};

/* Takes ownership of context, including on failure; borrows fapl and ops. */
h5c_status_t h5c__viz_open(const char* path, double time, hid_t fapl, const h5c_viz_ops* ops, void* context, h5c_viz_t** out);

#endif
