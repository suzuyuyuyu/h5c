/* Serial visualization writer coverage, including the documented layout. */
#include "h5c_test.h"

#include "h5c/h5c_viz.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

H5C_TEST_MAIN_STATE;

#define PATH     "test_viz.h5"
#define BAD_PATH "test_viz_bad.h5"
#define NPOINTS  5
#define NCELLS   2
#define NPE      3
#define NPARTS   4

static const double nodes[NPOINTS * 3] = {
    0.1, 10.2, 100.3, 1.1, 11.2, 101.3, 2.1, 12.2, 102.3,
    3.1, 13.2, 103.3, 4.1, 14.2, 104.3
};
static const int32_t conn[NCELLS * NPE] = { 0, 3, 1, 4, 2, 0 };

static void check_info(h5c_file_t *f, const char *path, int rank,
                       size_t d0, size_t d1, h5c_type_t type)
{
    h5c_dataset_info_t info;

    H5C_CHECK(h5c_dataset_info(f, path, &info));
    H5C_ASSERT(info.rank == rank, "%s rank: got %d, want %d", path,
               info.rank, rank);
    H5C_ASSERT_EQ_SIZE(info.dims[0], d0, path);
    if (rank == 2) {
        H5C_ASSERT_EQ_SIZE(info.dims[1], d1, path);
    }
    H5C_ASSERT(info.type == type, "%s type: got %d, want %d", path,
               (int)info.type, (int)type);
}

static void check_attr_string(h5c_file_t *f, const char *path,
                              const char *name, const char *want)
{
    char *got = NULL;

    H5C_CHECK(h5c_read_attr_str(f, path, name, &got));
    H5C_ASSERT(got != NULL && strcmp(got, want) == 0,
               "%s/%s: got '%s', want '%s'", path, name,
               got != NULL ? got : "(null)", want);
    h5c_free_string(got);
}

static void check_attribute_type(h5c_file_t *f, const char *path,
                                 const char *want)
{
    check_attr_string(f, path, "attribute_type", want);
}

static void write_fixture(h5c_viz_t *viz)
{
    h5c_viz_mesh_t mesh;
    double p[NPOINTS], stress[NPOINTS * 6];
    float velocity[3][NPOINTS], cell_velocity[3][NCELLS];
    int32_t subdomain[NCELLS];
    double cell_stress[NCELLS * 6];
    double px[NPARTS], py[NPARTS], pz[NPARTS], radius[NPARTS];
    const void *point_comps[3], *cell_comps[3], *xyz[3];
    size_t i, c;

    memset(&mesh, 0, sizeof mesh);
    mesh.kind = H5C_VIZ_UNSTRUCTURED;
    mesh.name = "fluid";
    mesh.topology = "Triangle";
    mesh.nodes_per_element = NPE;
    mesh.num_points = NPOINTS;
    mesh.num_cells = NCELLS;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));

    for (i = 0; i < NPOINTS; i++) {
        p[i] = 10.0 + 0.25 * (double)i;
        for (c = 0; c < 3; c++) {
            velocity[c][i] = (float)(100.0 * (double)(c + 1) + i + 0.5);
        }
        for (c = 0; c < 6; c++) {
            stress[i * 6 + c] = 1000.0 * (double)(c + 1) + i + 0.75;
        }
    }
    for (i = 0; i < NCELLS; i++) {
        subdomain[i] = (int32_t)(20 + i * 7);
        for (c = 0; c < 3; c++) {
            cell_velocity[c][i] = (float)(200.0 * (double)(c + 1) + i + 0.25);
        }
        for (c = 0; c < 6; c++) {
            cell_stress[i * 6 + c] = 3000.0 * (double)(c + 1) + i + 0.5;
        }
    }

    H5C_CHECK(h5c_viz_write_nodes(viz, nodes, H5C_F64));
    H5C_CHECK(h5c_viz_write_connectivity(viz, conn, H5C_I32));
    H5C_CHECK(h5c_viz_write_point_data(viz, "Pressure", p, H5C_F64, 1));
    point_comps[0] = velocity[0];
    point_comps[1] = velocity[1];
    point_comps[2] = velocity[2];
    H5C_CHECK(h5c_viz_write_point_data_comps(viz, "Velocity", point_comps,
                                             H5C_F32, 3));
    H5C_CHECK(h5c_viz_write_point_data(viz, "Stress", stress, H5C_F64, 6));
    H5C_CHECK(h5c_viz_write_cell_data(viz, "Subdomain", subdomain, H5C_I32, 1));
    cell_comps[0] = cell_velocity[0];
    cell_comps[1] = cell_velocity[1];
    cell_comps[2] = cell_velocity[2];
    H5C_CHECK(h5c_viz_write_cell_data_comps(viz, "CellVelocity", cell_comps,
                                            H5C_F32, 3));
    H5C_CHECK(h5c_viz_write_cell_data(viz, "CellStress", cell_stress,
                                     H5C_F64, 6));

    memset(&mesh, 0, sizeof mesh);
    mesh.kind = H5C_VIZ_POLYDATA;
    mesh.name = "particles";
    mesh.num_points = NPARTS;
    for (i = 0; i < NPARTS; i++) {
        px[i] = 0.5 + 2.0 * i;
        py[i] = 10.5 + 3.0 * i;
        pz[i] = 20.5 + 5.0 * i;
        radius[i] = 0.01 + 0.1 * i;
    }
    xyz[0] = px;
    xyz[1] = py;
    xyz[2] = pz;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));
    H5C_CHECK(h5c_viz_write_nodes_comps(viz, xyz, H5C_F64));
    H5C_CHECK(h5c_viz_write_point_data(viz, "Radius", radius, H5C_F64, 1));
}

static h5c_status_t open_serial_fixture(const char *path)
{
    h5c_viz_t *viz = NULL;
    h5c_status_t st = h5c_viz_open(path, 12.5, &viz);

    if (st == H5C_OK) {
        write_fixture(viz);
        st = h5c_viz_status(viz);
        if (h5c_viz_close(viz) != H5C_OK && st == H5C_OK) {
            st = H5C_ERR_HDF5;
        }
    }
    return st;
}

static void check_fixture(const char *path)
{
    h5c_file_t *f = NULL;
    int32_t scheme = -1, npe = -1;
    double time = -1.0;
    size_t poff = 99, coff = 99;
    double got_nodes[NPOINTS * 3], got_p[NPOINTS], got_stress[NPOINTS * 6];
    float got_velocity[NPOINTS * 3], got_cell_velocity[NCELLS * 3];
    int32_t got_conn[NCELLS * NPE], got_subdomain[NCELLS];
    double got_cell_stress[NCELLS * 6], got_radius[NPARTS];
    size_t dims_nodes[2] = { NPOINTS, 3 };
    size_t dims_conn[2] = { NCELLS, NPE };
    size_t dims_velocity[2] = { NPOINTS, 3 };
    size_t dims_stress[2] = { NPOINTS, 6 };
    size_t dims_cell_velocity[2] = { NCELLS, 3 };
    size_t dims_cell_stress[2] = { NCELLS, 6 };
    const size_t dims_particles[2] = { NPARTS, 3 };
    size_t i, c;

    H5C_CHECK(h5c_open(path, H5C_READ, &f));
    if (f == NULL) {
        return;
    }
    H5C_CHECK(h5c_read_attr_scalar(f, "/", "scheme_version", &scheme,
                                   H5C_I32));
    H5C_ASSERT(scheme == H5C_SCHEME_VERSION, "scheme_version: got %d", scheme);
    H5C_CHECK(h5c_read_attr_scalar(f, "/", "time", &time, H5C_F64));
    H5C_ASSERT(time == 12.5, "time: got %.17g", time);
    check_attr_string(f, "/fluid", "topology_type", "Triangle");
    H5C_CHECK(h5c_read_attr_scalar(f, "/fluid", "nodes_per_element", &npe,
                                   H5C_I32));
    H5C_ASSERT(npe == NPE, "nodes_per_element: got %d", npe);

    check_info(f, "/fluid/geometry/nodes", 2, NPOINTS, 3, H5C_F64);
    check_info(f, "/fluid/geometry/connectivity", 2, NCELLS, NPE, H5C_I32);
    check_info(f, "/fluid/point_data/Pressure", 1, NPOINTS, 0, H5C_F64);
    check_info(f, "/fluid/point_data/Velocity", 2, NPOINTS, 3, H5C_F32);
    check_info(f, "/fluid/point_data/Stress", 2, NPOINTS, 6, H5C_F64);
    check_info(f, "/fluid/cell_data/Subdomain", 1, NCELLS, 0, H5C_I32);
    check_info(f, "/fluid/cell_data/CellVelocity", 2, NCELLS, 3, H5C_F32);
    check_info(f, "/fluid/cell_data/CellStress", 2, NCELLS, 6, H5C_F64);
    check_attribute_type(f, "/fluid/point_data/Pressure", "Scalar");
    check_attribute_type(f, "/fluid/point_data/Velocity", "Vector");
    check_attribute_type(f, "/fluid/point_data/Stress", "Tensor6");
    check_attribute_type(f, "/fluid/cell_data/Subdomain", "Scalar");
    check_attribute_type(f, "/fluid/cell_data/CellVelocity", "Vector");
    check_attribute_type(f, "/fluid/cell_data/CellStress", "Tensor6");

    H5C_CHECK(h5c_read(f, "/fluid/geometry/nodes", got_nodes, H5C_F64,
                       2, dims_nodes));
    H5C_CHECK(h5c_read(f, "/fluid/geometry/connectivity", got_conn, H5C_I32,
                       2, dims_conn));
    H5C_CHECK(h5c_read(f, "/fluid/point_data/Pressure", got_p, H5C_F64,
                       1, &dims_nodes[0]));
    H5C_CHECK(h5c_read(f, "/fluid/point_data/Velocity", got_velocity,
                       H5C_F32, 2, dims_velocity));
    H5C_CHECK(h5c_read(f, "/fluid/point_data/Stress", got_stress, H5C_F64,
                       2, dims_stress));
    H5C_CHECK(h5c_read(f, "/fluid/cell_data/Subdomain", got_subdomain,
                       H5C_I32, 1, &dims_conn[0]));
    H5C_CHECK(h5c_read(f, "/fluid/cell_data/CellVelocity", got_cell_velocity,
                       H5C_F32, 2, dims_cell_velocity));
    H5C_CHECK(h5c_read(f, "/fluid/cell_data/CellStress", got_cell_stress,
                       H5C_F64, 2, dims_cell_stress));
    /* The component-array API stores rows as p0c0 p0c1 p0c2 p1c0 ... . */
    for (i = 0; i < NPOINTS; i++) {
        for (c = 0; c < 3; c++) {
            H5C_ASSERT(got_velocity[i * 3 + c] ==
                           (float)(100.0 * (c + 1) + i + 0.5),
                       "raw velocity[%lu * 3 + %lu] mismatch",
                       (unsigned long)i, (unsigned long)c);
        }
    }
    for (i = 0; i < NCELLS; i++) {
        for (c = 0; c < 3; c++) {
            H5C_ASSERT(got_cell_velocity[i * 3 + c] ==
                           (float)(200.0 * (c + 1) + i + 0.25),
                       "raw cell velocity[%lu * 3 + %lu] mismatch",
                       (unsigned long)i, (unsigned long)c);
        }
    }
    for (i = 0; i < NPOINTS * 3; i++) {
        H5C_ASSERT(got_nodes[i] == nodes[i], "nodes[%lu] mismatch",
                   (unsigned long)i);
    }
    for (i = 0; i < NCELLS * NPE; i++) {
        H5C_ASSERT(got_conn[i] == conn[i], "connectivity[%lu] mismatch",
                   (unsigned long)i);
    }
    for (i = 0; i < NPOINTS; i++) {
        H5C_ASSERT(got_p[i] == 10.0 + 0.25 * i, "pressure[%lu] mismatch",
                   (unsigned long)i);
        for (c = 0; c < 3; c++) {
            H5C_ASSERT(got_velocity[i * 3 + c] ==
                           (float)(100.0 * (c + 1) + i + 0.5),
                       "velocity[%lu][%lu] mismatch", (unsigned long)i,
                       (unsigned long)c);
        }
        for (c = 0; c < 6; c++) {
            H5C_ASSERT(got_stress[i * 6 + c] == 1000.0 * (c + 1) + i + 0.75,
                       "stress[%lu][%lu] mismatch", (unsigned long)i,
                       (unsigned long)c);
        }
    }
    for (i = 0; i < NCELLS; i++) {
        H5C_ASSERT(got_subdomain[i] == (int32_t)(20 + i * 7),
                   "subdomain[%lu] mismatch", (unsigned long)i);
        for (c = 0; c < 3; c++) {
            H5C_ASSERT(got_cell_velocity[i * 3 + c] ==
                           (float)(200.0 * (c + 1) + i + 0.25),
                       "cell velocity[%lu][%lu] mismatch", (unsigned long)i,
                       (unsigned long)c);
        }
        for (c = 0; c < 6; c++) {
            H5C_ASSERT(got_cell_stress[i * 6 + c] ==
                           3000.0 * (c + 1) + i + 0.5,
                       "cell stress[%lu][%lu] mismatch", (unsigned long)i,
                       (unsigned long)c);
        }
    }

    check_attr_string(f, "/particles", "topology_type", "Polyvertex");
    H5C_CHECK(h5c_read_attr_scalar(f, "/particles", "nodes_per_element",
                                   &npe, H5C_I32));
    H5C_ASSERT(npe == 1, "polydata nodes_per_element: got %d", npe);
    H5C_ASSERT(!h5c_exists(f, "/particles/geometry/connectivity"),
               "polydata unexpectedly has connectivity");
    H5C_ASSERT(!h5c_exists(f, "/particles/cell_data"),
               "polydata unexpectedly has cell_data");
    check_info(f, "/particles/geometry/nodes", 2, NPARTS, 3, H5C_F64);
    check_info(f, "/particles/point_data/Radius", 1, NPARTS, 0, H5C_F64);
    H5C_CHECK(h5c_read(f, "/particles/point_data/Radius", got_radius,
                       H5C_F64, 1, &dims_particles[0]));
    for (i = 0; i < NPARTS; i++) {
        H5C_ASSERT(got_radius[i] == 0.01 + 0.1 * i,
                   "radius[%lu] mismatch", (unsigned long)i);
    }
    H5C_CHECK(h5c_close(f));
    H5C_ASSERT_EQ_SIZE(poff, 99, "unused serial point offset sentinel");
    H5C_ASSERT_EQ_SIZE(coff, 99, "unused serial cell offset sentinel");
}

static void check_offsets(void)
{
    h5c_viz_t *viz = NULL;
    h5c_viz_mesh_t mesh = {0};
    size_t poff = 99, coff = 99;

    H5C_CHECK(h5c_viz_open("test_viz_offsets.h5", 0.0, &viz));
    if (viz == NULL) {
        return;
    }
    mesh.kind = H5C_VIZ_UNSTRUCTURED;
    mesh.num_points = NPOINTS;
    mesh.num_cells = NCELLS;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));
    H5C_CHECK(h5c_viz_offsets(viz, &poff, &coff));
    H5C_ASSERT_EQ_SIZE(poff, 0, "serial point offset");
    H5C_ASSERT_EQ_SIZE(coff, 0, "serial cell offset");
    H5C_CHECK(h5c_viz_close(viz));
}

static void check_connectivity_type(const char *path, h5c_type_t type)
{
    h5c_viz_t *viz = NULL;
    h5c_file_t *f = NULL;
    h5c_viz_mesh_t mesh = {0};
    const size_t dims[2] = { NCELLS, NPE };
    int8_t c8[NCELLS * NPE];
    int16_t c16[NCELLS * NPE];
    int32_t c32[NCELLS * NPE];
    int64_t c64[NCELLS * NPE];
    size_t i;

    H5C_CHECK(h5c_viz_open(path, 0.0, &viz));
    if (viz == NULL) {
        return;
    }
    mesh.kind = H5C_VIZ_UNSTRUCTURED;
    mesh.name = "mesh";
    mesh.nodes_per_element = NPE;
    mesh.num_points = NPOINTS;
    mesh.num_cells = NCELLS;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));
    switch (type) {
    case H5C_I8:
        for (i = 0; i < NCELLS * NPE; i++) c8[i] = (int8_t)conn[i];
        H5C_CHECK(h5c_viz_write_connectivity(viz, c8, type));
        break;
    case H5C_I16:
        for (i = 0; i < NCELLS * NPE; i++) c16[i] = (int16_t)conn[i];
        H5C_CHECK(h5c_viz_write_connectivity(viz, c16, type));
        break;
    case H5C_I32:
        H5C_CHECK(h5c_viz_write_connectivity(viz, conn, type));
        break;
    default:
        for (i = 0; i < NCELLS * NPE; i++) c64[i] = conn[i];
        H5C_CHECK(h5c_viz_write_connectivity(viz, c64, type));
        break;
    }
    H5C_CHECK(h5c_viz_close(viz));
    H5C_CHECK(h5c_open(path, H5C_READ, &f));
    if (f == NULL) {
        return;
    }
    check_info(f, "/mesh/geometry/connectivity", 2, NCELLS, NPE, type);
    switch (type) {
    case H5C_I8:  H5C_CHECK(h5c_read(f, "/mesh/geometry/connectivity", c8,
                                    type, 2, dims)); break;
    case H5C_I16: H5C_CHECK(h5c_read(f, "/mesh/geometry/connectivity", c16,
                                    type, 2, dims)); break;
    case H5C_I32: H5C_CHECK(h5c_read(f, "/mesh/geometry/connectivity", c32,
                                    type, 2, dims)); break;
    default:      H5C_CHECK(h5c_read(f, "/mesh/geometry/connectivity", c64,
                                    type, 2, dims)); break;
    }
    for (i = 0; i < NCELLS * NPE; i++) {
        int64_t got;

        switch (type) {
        case H5C_I8:  got = c8[i];  break;
        case H5C_I16: got = c16[i]; break;
        case H5C_I32: got = c32[i]; break;
        default:      got = c64[i]; break;
        }
        H5C_ASSERT(got == conn[i], "connectivity type %d[%lu] mismatch",
                   (int)type, (unsigned long)i);
    }
    H5C_CHECK(h5c_close(f));
}

static void check_bad_connectivity(void)
{
    h5c_viz_t *viz = NULL;
    h5c_file_t *f = NULL;
    h5c_viz_mesh_t mesh = {0};
    int32_t bad[NCELLS * NPE] = { 0, 1, NPOINTS, 3, 4, 0 };
    h5c_status_t st;

    H5C_CHECK(h5c_viz_open(BAD_PATH, 0.0, &viz));
    if (viz == NULL) {
        return;
    }
    mesh.kind = H5C_VIZ_UNSTRUCTURED;
    mesh.name = "bad";
    mesh.nodes_per_element = NPE;
    mesh.num_points = NPOINTS;
    mesh.num_cells = NCELLS;
    H5C_CHECK(h5c_viz_begin_mesh(viz, &mesh));
    st = h5c_viz_write_connectivity(viz, bad, H5C_I32);
    H5C_ASSERT(st == H5C_ERR_INVALID_ARG, "bad connectivity: got %s",
               h5c_status_string(st));
    H5C_CHECK(h5c_viz_close(viz));
    H5C_CHECK(h5c_open(BAD_PATH, H5C_READ, &f));
    if (f != NULL) {
        H5C_ASSERT(!h5c_exists(f, "/bad/geometry/connectivity"),
                   "rejected connectivity was written");
        H5C_CHECK(h5c_close(f));
    }
}

int main(void)
{
    H5C_CHECK(h5c_init());
    H5C_CHECK(open_serial_fixture(PATH));
    check_fixture(PATH);
    check_offsets();
    check_connectivity_type("test_viz_i8.h5", H5C_I8);
    check_connectivity_type("test_viz_i16.h5", H5C_I16);
    check_connectivity_type("test_viz_i32.h5", H5C_I32);
    check_connectivity_type("test_viz_i64.h5", H5C_I64);
    check_bad_connectivity();

    h5c_finalize();
    return H5C_TEST_SUMMARY("test_viz");
}
