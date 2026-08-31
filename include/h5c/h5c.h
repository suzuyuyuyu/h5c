/*
 * h5c — a compact C wrapper around the HDF5 C API.
 *
 * Design contract (see docs/FORMAT.md and docs/USAGE.md):
 *
 *   - Arrays are FLAT buffers plus an explicit `dims` array.
 *   - Dimension order is ROW-MAJOR: dims[rank-1] varies fastest.
 *     This matches the HDF5 C API, so no index translation happens here.
 *     A Fortran array a(nx, ny) corresponds to dims = {ny, nx}.
 *   - Every function returns h5c_status_t. H5C_OK (0) means success.
 *   - The file handle also keeps a STICKY status: the first non-OK status
 *     seen on that file. A later success never clears it.
 *
 * This header never includes <mpi.h>. Parallel I/O lives in h5c_mpi.h.
 */
#ifndef H5C_H
#define H5C_H

#include <stddef.h>
#include <stdint.h>

#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum rank handled by the info struct. HDF5's own limit is 32. */
#define H5C_MAX_RANK 32

/* ------------------------------------------------------------------ */
/* Status                                                              */
/* ------------------------------------------------------------------ */

/*
 * Numeric values are APPEND-ONLY. Never renumber an existing member:
 * callers may have been compiled against an older header.
 */
typedef enum h5c_status {
    H5C_OK              = 0,
    H5C_ERR_INVALID_ARG = 1,  /* NULL pointer, bad rank, unknown enum value  */
    H5C_ERR_NOT_FOUND   = 2,  /* no such file, group, dataset or attribute   */
    H5C_ERR_SHAPE_MISMATCH = 3, /* rank or extents differ from the dataset   */
    H5C_ERR_TYPE_MISMATCH  = 4, /* datatype class cannot be converted        */
    H5C_ERR_EXISTS      = 5,  /* object exists and replace was not requested */
    H5C_ERR_HDF5        = 6,  /* the HDF5 library reported a failure         */
    H5C_ERR_MPI         = 7,  /* an MPI call failed                          */
    H5C_ERR_NOMEM       = 8,  /* allocation failed                           */
    H5C_ERR_STATE       = 9,  /* wrong lifecycle state (e.g. read-only file) */
    H5C_ERR_UNSUPPORTED = 10  /* valid request this build cannot serve       */
} h5c_status_t;

/* Static, never-NULL description of a status value. */
const char *h5c_status_string(h5c_status_t status);

/*
 * Details of the most recent failure on THIS THREAD. The returned pointer is
 * owned by h5c and stays valid until the next failing call on the same thread.
 * `hdf5_err` is the raw herr_t when status is H5C_ERR_HDF5, otherwise 0.
 */
typedef struct h5c_error {
    h5c_status_t status;
    long         hdf5_err;
    char         message[256];
} h5c_error_t;

const h5c_error_t *h5c_last_error(void);

/* ------------------------------------------------------------------ */
/* Library lifecycle                                                   */
/* ------------------------------------------------------------------ */

/*
 * Optional. Every entry point initialises the library on first use, so callers
 * may skip this. Its real job is installing h5c's error handling, which
 * suppresses HDF5's own multi-page stack dump on stderr.
 */
h5c_status_t h5c_init(void);

/* Releases h5c's cached datatypes and calls H5close(). Optional. */
void h5c_finalize(void);

/*
 * 0 = quiet, h5c reports errors only through the return value and
 *     h5c_last_error() (the default).
 * 1 = also let HDF5 print its native error stack to stderr.
 */
void h5c_set_error_verbosity(int level);

/* ------------------------------------------------------------------ */
/* Datatypes                                                           */
/* ------------------------------------------------------------------ */

typedef enum h5c_type {
    H5C_TYPE_UNKNOWN = 0,
    H5C_F32    = 1,
    H5C_F64    = 2,
    H5C_I32    = 3,
    H5C_I64    = 4,
    H5C_BOOL   = 5,
    H5C_STRING = 6
} h5c_type_t;

/*
 * Booleans are stored as an enum over int8 with members FALSE=0 and TRUE=1.
 * That is one byte per element, h5dump prints TRUE/FALSE, and h5py reads the
 * dataset back as numpy bool. h5fortran reads it as `logical` because its
 * logical reader asks for H5T_NATIVE_INTEGER and HDF5 converts.
 *
 * Values other than 0 and 1 are written verbatim and are NOT validated:
 * normalising would cost a full extra pass over the buffer.
 */
typedef int8_t h5c_bool_t;
#define H5C_FALSE ((h5c_bool_t)0)
#define H5C_TRUE  ((h5c_bool_t)1)

/* Size of one element of `type` in memory. Returns 0 for unknown/string. */
size_t h5c_type_size(h5c_type_t type);

/* ------------------------------------------------------------------ */
/* Files                                                               */
/* ------------------------------------------------------------------ */

typedef struct h5c_file h5c_file_t;

typedef enum h5c_mode {
    H5C_READ      = 1,  /* open an existing file read-only        */
    H5C_READWRITE = 2,  /* open an existing file for reading and writing */
    H5C_TRUNCATE  = 3   /* create, discarding any existing file   */
} h5c_mode_t;

h5c_status_t h5c_open(const char *path, h5c_mode_t mode, h5c_file_t **out);

/*
 * Closes the file and frees the handle, even on failure: the handle is invalid
 * after this call regardless of the returned status.
 *
 * The returned status describes ONLY the close itself, so that a caller can
 * distinguish "the file may not have been flushed" from "an earlier call
 * failed and was handled". The sticky status is NOT reported here; read it
 * with h5c_file_status() before closing:
 *
 *     h5c_status_t failed = h5c_file_status(f);  // did anything go wrong?
 *     h5c_status_t closed = h5c_close(f);        // did the close succeed?
 */
h5c_status_t h5c_close(h5c_file_t *file);

/* First non-OK status recorded on this file, or H5C_OK. */
h5c_status_t h5c_file_status(const h5c_file_t *file);
void         h5c_file_clear_status(h5c_file_t *file);

/* Interop: the underlying HDF5 id, or H5I_INVALID_HID. Borrowed, do not close. */
hid_t h5c_file_hid(const h5c_file_t *file);

/*
 * Wraps an hid_t that the caller continues to own. h5c_close() on the result
 * frees the wrapper but does NOT close `fid`.
 */
h5c_status_t h5c_file_from_hid(hid_t fid, h5c_file_t **out);

/*
 * True (1) when `file` was opened through the parallel entry points in
 * h5c_mpi.h, so that generic code can branch on it WITHOUT including <mpi.h>.
 * Always 0 in a build without parallel support, and 0 for a NULL handle.
 */
int h5c_is_parallel(const h5c_file_t *file);

/* ------------------------------------------------------------------ */
/* Datasets                                                            */
/* ------------------------------------------------------------------ */

/* Flags for h5c_write(). */
#define H5C_WRITE_DEFAULT 0u
#define H5C_WRITE_REPLACE 1u  /* delete an existing dataset first */

typedef struct h5c_dataset_info {
    int        rank;
    size_t     dims[H5C_MAX_RANK];
    h5c_type_t type;   /* H5C_TYPE_UNKNOWN if the stored type has no mapping */
    size_t     count;  /* product of dims; 1 for a scalar dataset */
} h5c_dataset_info_t;

/*
 * Writes `buf` to `path`, creating intermediate groups as needed.
 *
 * rank == 0 writes a scalar dataset; `dims` may then be NULL.
 * Without H5C_WRITE_REPLACE an existing dataset is written in place, which
 * requires its rank, extents and type to match; otherwise H5C_ERR_EXISTS or
 * H5C_ERR_SHAPE_MISMATCH is returned.
 *
 * A zero extent is legal: it creates an empty dataset of that shape and
 * transfers nothing, and `buf` may then be NULL. Parallel I/O requires this
 * (a rank may own no rows), and the same rule applies in serial so that the
 * library has exactly one validation policy.
 */
h5c_status_t h5c_write(h5c_file_t *file, const char *path, const void *buf,
                       h5c_type_t type, int rank, const size_t *dims,
                       unsigned flags);

/*
 * Reads `path` into `buf`, which the caller has already sized.
 * The dataset's rank and extents must equal `rank` and `dims` exactly.
 * HDF5 converts between the stored and requested types where it can.
 *
 * As for h5c_write(), a zero extent is legal and transfers nothing.
 */
h5c_status_t h5c_read(h5c_file_t *file, const char *path, void *buf,
                      h5c_type_t type, int rank, const size_t *dims);

h5c_status_t h5c_dataset_info(h5c_file_t *file, const char *path,
                              h5c_dataset_info_t *out);

/*
 * True (1) when `path` resolves to an existing object in `file`.
 *
 * This is a pure query: it never touches the file's sticky status, and it
 * returns 0 both for "absent" and for "the handle or path was unusable".
 * Callers that must distinguish the two should ask h5c_dataset_info() and
 * inspect the returned status instead.
 */
int h5c_exists(h5c_file_t *file, const char *path);

/*
 * Convenience for pre- and post-processing: allocates a buffer of the right
 * size, reads into it, and reports the shape. Release it with h5c_free().
 * Solver code should prefer h5c_dataset_info() plus h5c_read().
 */
h5c_status_t h5c_read_alloc(h5c_file_t *file, const char *path,
                            h5c_type_t type, void **buf,
                            h5c_dataset_info_t *info);

/*
 * For an empty dataset h5c_read_alloc() still yields a non-NULL, freeable
 * pointer, so that `*buf == NULL` unambiguously means failure.
 */

void h5c_free(void *ptr);

/* ------------------------------------------------------------------ */
/* Strings                                                             */
/* ------------------------------------------------------------------ */

/*
 * Strings are stored as a scalar dataset of fixed-length H5T_C_S1, which is
 * what h5fortran writes and reads. The stored length is strlen(value).
 */
h5c_status_t h5c_write_string(h5c_file_t *file, const char *path,
                              const char *value, unsigned flags);

/*
 * Variable-length alternative. Convenient from C, but h5fortran cannot read
 * it: its character reader expects a fixed-length type.
 */
h5c_status_t h5c_write_string_vlen(h5c_file_t *file, const char *path,
                                   const char *value, unsigned flags);

/*
 * Reads either representation. `*out` is allocated by h5c and is always
 * NUL-terminated; release it with h5c_free_string().
 */
h5c_status_t h5c_read_string(h5c_file_t *file, const char *path, char **out);

void h5c_free_string(char *s);

/* ------------------------------------------------------------------ */
/* Attributes                                                          */
/* ------------------------------------------------------------------ */

/*
 * `obj_path` may name a dataset, a group, or "/" for the root group.
 * Attributes are written separately from the dataset they annotate; there is
 * no combined "write data and attributes" call.
 */
h5c_status_t h5c_write_attr_str(h5c_file_t *file, const char *obj_path,
                                const char *name, const char *value);

/* `*out` is allocated by h5c; release it with h5c_free_string(). */
h5c_status_t h5c_read_attr_str(h5c_file_t *file, const char *obj_path,
                               const char *name, char **out);

/* Scalar numeric attributes, e.g. "time", "step", "dt". */
h5c_status_t h5c_write_attr_scalar(h5c_file_t *file, const char *obj_path,
                                   const char *name, const void *value,
                                   h5c_type_t type);

h5c_status_t h5c_read_attr_scalar(h5c_file_t *file, const char *obj_path,
                                  const char *name, void *value,
                                  h5c_type_t type);

int h5c_attr_exists(h5c_file_t *file, const char *obj_path, const char *name);

/* ------------------------------------------------------------------ */
/* Interleaved multi-component fields                                  */
/* ------------------------------------------------------------------ */

/*
 * Vector and tensor fields are stored INTERLEAVED as [n, ncomp], because that
 * is what XDMF3 Vector/Tensor attributes require and what h5fortran's
 * visualization writer produces from field(ncomp, nlocal).
 *
 * Solvers commonly hold components separately (u[n], v[n], w[n]). These calls
 * take an array of `ncomp` component pointers and handle the interleaving.
 *
 * ncomp is 1, 3, 6 or 9 for Scalar, Vector, Tensor6 and Tensor respectively.
 * For ncomp == 6 the XDMF component order is XX, XY, XZ, YY, YZ, ZZ.
 *
 * WRITING PACKS, READING STRIDES. Writing gathers the components into a
 * contiguous staging buffer and issues one write, because a strided
 * collective write forces read-modify-write in MPI-IO. Reading one component
 * has no such penalty, so h5c_read_component() selects a strided hyperslab
 * and moves only the bytes that component needs.
 */
h5c_status_t h5c_write_interleaved(h5c_file_t *file, const char *path,
                                   const void *const *comps, size_t ncomp,
                                   size_t n, h5c_type_t type, unsigned flags);

h5c_status_t h5c_read_interleaved(h5c_file_t *file, const char *path,
                                  void *const *comps, size_t ncomp,
                                  size_t n, h5c_type_t type);

/* Reads component `comp` (0-based) of an [n, ncomp] dataset into `buf`. */
h5c_status_t h5c_read_component(h5c_file_t *file, const char *path,
                                void *buf, size_t comp, size_t n,
                                h5c_type_t type);

/*
 * Upper bound on the staging buffer used when packing, in bytes. Larger
 * fields are written in row-wise tiles, so memory stays bounded while every
 * write remains contiguous in the file. Default 256 MiB.
 */
void   h5c_set_pack_limit(size_t bytes);
size_t h5c_pack_limit(void);

/* ------------------------------------------------------------------ */
/* Typed convenience wrappers                                          */
/* ------------------------------------------------------------------ */

/*
 * One set per type: scalar, 1-D, and N-D. Rank is carried by `dims`, so there
 * is no per-rank expansion. Use h5c_write() directly when you need flags.
 */
#define H5C_DECLARE_TYPED_IO(SUF, CTYPE, ENUMV)                               \
    static inline h5c_status_t h5c_write_##SUF(                               \
        h5c_file_t *f, const char *p, const CTYPE *buf,                       \
        int rank, const size_t *dims) {                                       \
        return h5c_write(f, p, buf, ENUMV, rank, dims, H5C_WRITE_DEFAULT);    \
    }                                                                         \
    static inline h5c_status_t h5c_write_##SUF##_1d(                          \
        h5c_file_t *f, const char *p, const CTYPE *buf, size_t n) {           \
        return h5c_write(f, p, buf, ENUMV, 1, &n, H5C_WRITE_DEFAULT);         \
    }                                                                         \
    static inline h5c_status_t h5c_write_##SUF##_scalar(                      \
        h5c_file_t *f, const char *p, CTYPE value) {                          \
        return h5c_write(f, p, &value, ENUMV, 0, NULL, H5C_WRITE_DEFAULT);    \
    }                                                                         \
    static inline h5c_status_t h5c_read_##SUF(                                \
        h5c_file_t *f, const char *p, CTYPE *buf,                             \
        int rank, const size_t *dims) {                                       \
        return h5c_read(f, p, buf, ENUMV, rank, dims);                        \
    }                                                                         \
    static inline h5c_status_t h5c_read_##SUF##_1d(                           \
        h5c_file_t *f, const char *p, CTYPE *buf, size_t n) {                 \
        return h5c_read(f, p, buf, ENUMV, 1, &n);                             \
    }                                                                         \
    static inline h5c_status_t h5c_read_##SUF##_scalar(                       \
        h5c_file_t *f, const char *p, CTYPE *value) {                         \
        return h5c_read(f, p, value, ENUMV, 0, NULL);                         \
    }                                                                         \
    /* Interleaved forms. These also spare the caller the cast to             \
       `const void *const *`, which C requires but which reads badly. */      \
    static inline h5c_status_t h5c_write_interleaved_##SUF(                   \
        h5c_file_t *f, const char *p, const CTYPE *const *comps,              \
        size_t ncomp, size_t n) {                                             \
        return h5c_write_interleaved(f, p, (const void *const *)comps,        \
                                     ncomp, n, ENUMV, H5C_WRITE_DEFAULT);     \
    }                                                                         \
    static inline h5c_status_t h5c_read_interleaved_##SUF(                    \
        h5c_file_t *f, const char *p, CTYPE *const *comps,                    \
        size_t ncomp, size_t n) {                                             \
        return h5c_read_interleaved(f, p, (void *const *)comps,               \
                                    ncomp, n, ENUMV);                         \
    }                                                                         \
    static inline h5c_status_t h5c_read_component_##SUF(                      \
        h5c_file_t *f, const char *p, CTYPE *buf, size_t comp, size_t n) {    \
        return h5c_read_component(f, p, buf, comp, n, ENUMV);                 \
    }

H5C_DECLARE_TYPED_IO(f32,  float,      H5C_F32)
H5C_DECLARE_TYPED_IO(f64,  double,     H5C_F64)
H5C_DECLARE_TYPED_IO(i32,  int32_t,    H5C_I32)
H5C_DECLARE_TYPED_IO(i64,  int64_t,    H5C_I64)
H5C_DECLARE_TYPED_IO(bool, h5c_bool_t, H5C_BOOL)

#ifdef __cplusplus
}
#endif

#endif /* H5C_H */
