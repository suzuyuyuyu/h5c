/* Internal helpers shared by the h5c translation units. Not installed. */
#ifndef H5C_INTERNAL_H
#define H5C_INTERNAL_H

#include "h5c/h5c.h"

#ifdef H5C_HAVE_PARALLEL
#  include <mpi.h>
#endif

#include <stdio.h>

#if defined(__cplusplus)
#  define H5C_THREAD_LOCAL thread_local
#elif defined(_MSC_VER)
#  define H5C_THREAD_LOCAL __declspec(thread)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_THREADS__)
#  define H5C_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
#  define H5C_THREAD_LOCAL __thread
#else
#  define H5C_THREAD_LOCAL /* single-threaded fallback */
#endif

struct h5c_file {
    hid_t        fid;
    h5c_status_t sticky;   /* first non-OK status seen on this file */
    int          borrowed; /* 1 when fid is owned by the caller */
    int          readonly;
#ifdef H5C_HAVE_PARALLEL
    int          parallel;   /* 1 when opened through the h5c_mpi.h entry points */
    int          collective; /* transfer mode; 1 (collective) by default */
    MPI_Comm     comm;       /* borrowed from the caller, valid until close */
#endif
};

/* --- error plumbing ------------------------------------------------ */

/* Records `status` with a printf-style message and returns `status`. */
h5c_status_t h5c__fail(h5c_status_t status, const char *fmt, ...);

/* As above, but also stores the raw HDF5 error code. */
h5c_status_t h5c__fail_hdf5(long herr, const char *fmt, ...);

/* Folds `status` into the file's sticky error and returns it unchanged. */
h5c_status_t h5c__record(h5c_file_t *file, h5c_status_t status);

/* Ensures the library is initialised. Safe to call repeatedly. */
h5c_status_t h5c__ensure_init(void);

/* --- datatype mapping ---------------------------------------------- */

/*
 * File and memory datatypes for `type`. Both are BORROWED: they belong to h5c
 * and must not be closed by the caller. Returns H5I_INVALID_HID on failure.
 *
 * For H5C_BOOL both are the same enum id, so no conversion path runs.
 */
hid_t h5c__file_type(h5c_type_t type);
hid_t h5c__mem_type(h5c_type_t type);

/*
 * Memory type to use when READING. Identical to h5c__mem_type() except for
 * H5C_BOOL, where it is H5T_NATIVE_INT8 rather than the boolean enum.
 *
 * HDF5 converts an enum to an integer but NOT an integer to an enum, so a
 * plain int dataset -- which is exactly what h5fortran writes for `logical` --
 * cannot be read through the enum. Reading through int8 works for both
 * representations, and h5c_bool_t is int8_t, so nothing else changes.
 *
 * Writing still uses the enum for both memory and file: that pairing needs no
 * conversion at all, and int8 -> enum is the direction HDF5 refuses.
 */
hid_t h5c__mem_type_read(h5c_type_t type);

/* Best-effort reverse mapping of a stored datatype. */
h5c_type_t h5c__type_from_hid(hid_t tid);

/* Releases the cached datatypes. Called by h5c_finalize(). */
void h5c__type_cleanup(void);

/* --- shared helpers ------------------------------------------------ */

/* Link creation property list that creates intermediate groups. */
hid_t h5c__lcpl(void);

/* --- dataset helpers, shared with the other translation units ------ */

/* Validates a file handle, path, rank and dims. */
h5c_status_t h5c__check_common(h5c_file_t *file, const char *path,
                               int rank, const size_t *dims);

/* Dataspace for rank/dims; rank 0 yields a scalar space. Caller closes it. */
hid_t h5c__make_space(int rank, const size_t *dims);

/* Product of dims; 1 for rank 0. Zero when any extent is zero. */
size_t h5c__count(int rank, const size_t *dims);

/* Fills `out` from an already-open dataset id. */
h5c_status_t h5c__info_from_dset(hid_t did, h5c_dataset_info_t *out);

/* Compares a stored shape against rank/dims, reporting a useful message. */
h5c_status_t h5c__shape_equals(const h5c_dataset_info_t *info,
                               const char *path, int rank, const size_t *dims);

/* --- interleave buffer arithmetic ----------------------------------- */

/*
 * The buffer-level half of the interleaved API, implemented once in
 * h5c_interleaved.c and used by both the serial and the parallel entry
 * points. Only the arithmetic is shared: the two public layers keep their own
 * validation, dataset handling and (in parallel) collective agreement.
 *
 * `esize` is the element size, `row0` the first row of the tile in the
 * caller's component arrays and `rows` the number of rows in it.
 */

/* Rows per tile: as many whole rows as h5c_pack_limit() allows, at least 1. */
size_t h5c__tile_rows(size_t n, size_t row_bytes);

/* Gathers rows [row0, row0+rows) of every component into `dst`. */
void h5c__pack_tile(char *dst, const void *const *comps, size_t ncomp,
                    size_t row0, size_t rows, size_t esize);

/* Scatters rows [row0, row0+rows) of `src` back into the components. */
void h5c__unpack_tile(const char *src, void *const *comps, size_t ncomp,
                      size_t row0, size_t rows, size_t esize);

#endif /* H5C_INTERNAL_H */
