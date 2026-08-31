/*
 * h5c — parallel (MPI) I/O.
 *
 * This header is the ONLY place h5c exposes <mpi.h>. Serial users include
 * h5c/h5c.h and never pull in MPI. It exists only when the library was built
 * with H5C_ENABLE_PARALLEL=ON, which also defines H5C_HAVE_PARALLEL.
 *
 * ---------------------------------------------------------------------------
 * Storage layout (identical to h5fortran, see docs/FORMAT.md)
 * ---------------------------------------------------------------------------
 *
 * Writing a distributed array to path P produces a GROUP:
 *
 *     P/data             all ranks' data concatenated along the split axis
 *     P/__partition__    int64 rank boundaries, length nprocs + 1
 *
 *     __partition__[0]      == 0
 *     rank r starts at         __partition__[r]
 *     rank r owns              __partition__[r+1] - __partition__[r] rows
 *     total rows            == __partition__[nprocs]
 *
 * The array is split along the SLOWEST-VARYING axis, dims[0] in C. That is the
 * same file axis h5fortran splits (its Fortran LAST dimension), so files
 * interoperate directly. A local extent of 0 is allowed. Every other
 * dimension must agree across all ranks.
 *
 * ---------------------------------------------------------------------------
 * Collective discipline
 * ---------------------------------------------------------------------------
 *
 * Every call here is collective over the file's communicator: all ranks must
 * call it, in the same order, with the same path. Transfers are collective by
 * default and only change when h5c_pset_collective() is called explicitly.
 *
 * Validation failures are agreed upon across ranks before any HDF5 collective
 * call, so a bad argument on one rank makes every rank return the same error
 * instead of deadlocking.
 */
#ifndef H5C_MPI_H
#define H5C_MPI_H

#include <mpi.h>

#include "h5c/h5c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Name of the rank-boundary dataset written beside "data". */
#define H5C_PARTITION_NAME "__partition__"

/* Opens `path` for parallel I/O on MPI_COMM_WORLD with MPI_INFO_NULL. */
h5c_status_t h5c_popen(const char *path, h5c_mode_t mode, h5c_file_t **out);

/*
 * As above on an explicit communicator. `info` carries MPI-IO hints and may
 * be MPI_INFO_NULL. The communicator must stay valid until h5c_close().
 */
h5c_status_t h5c_popen_comm(const char *path, h5c_mode_t mode,
                            MPI_Comm comm, MPI_Info info, h5c_file_t **out);

/*
 * Selects the transfer mode for subsequent parallel calls on this file.
 * Non-zero (the default) is collective. h5c never switches modes implicitly.
 *
 * This is NOT a collective call: it communicates with nobody and only mutates
 * per-file state. It is, however, state that must AGREE across ranks before
 * the next transfer, because a collective HDF5 transfer requires every rank to
 * have requested the same mode. Set it to the same value on every rank, or
 * ranks will disagree about a call they all have to enter together.
 */
h5c_status_t h5c_pset_collective(h5c_file_t *file, int collective);
int          h5c_pis_collective(const h5c_file_t *file);

/* The file's communicator, or MPI_COMM_NULL if it was not opened in parallel. */
MPI_Comm h5c_pcomm(const h5c_file_t *file);

/*
 * Writes this rank's block. `dims` describes the LOCAL block; dims[0] is this
 * rank's extent along the split axis and may be 0. dims[1..] must match on
 * every rank. rank >= 1 is required: a scalar has no axis to split.
 *
 * When this rank's block is empty (dims[0] == 0) `buf` may be NULL: there is
 * nothing to transfer, and the rank contributes an empty selection while still
 * entering the collective call. That is the normal shape of the code on a rank
 * that owns no rows, so it is guaranteed rather than merely tolerated.
 */
h5c_status_t h5c_pwrite(h5c_file_t *file, const char *path, const void *buf,
                        h5c_type_t type, int rank, const size_t *dims,
                        unsigned flags);

/*
 * Reads this rank's block, sized by `dims` exactly as for h5c_pwrite().
 * The stored __partition__ must have length nprocs + 1, start at 0, and be
 * non-decreasing, and its final value must equal the extent of "data" along
 * the split axis; otherwise the read fails before any data is transferred.
 *
 * As for h5c_pwrite(), `buf` may be NULL when dims[0] == 0.
 */
h5c_status_t h5c_pread(h5c_file_t *file, const char *path, void *buf,
                       h5c_type_t type, int rank, const size_t *dims);

/*
 * Shape of a distributed dataset. `local` reports this rank's block as
 * recorded in __partition__, `global` the full extent of "data".
 * Either output may be NULL.
 */
h5c_status_t h5c_pdataset_info(h5c_file_t *file, const char *path,
                               h5c_dataset_info_t *local,
                               h5c_dataset_info_t *global);

/*
 * Where this rank's block sits inside the distributed dataset.
 *
 * `*offset` is the start index along the split axis and `*nlocal` the number
 * of rows this rank owns; either output may be NULL. This is the one thing
 * h5c_pdataset_info() cannot tell you, and it exists so that callers never
 * have to open "<path>/" H5C_PARTITION_NAME and interpret the layout
 * themselves.
 *
 * Collective: every rank must call it, and the stored __partition__ is
 * validated exactly as h5c_pread() validates it.
 */
h5c_status_t h5c_poffset(h5c_file_t *file, const char *path,
                         size_t *offset, size_t *nlocal);

/*
 * The full rank-boundary vector: nprocs + 1 entries, starting at 0 and
 * non-decreasing (see docs/FORMAT.md).
 *
 * Pass `bounds == NULL` to learn the required length first, in which case
 * `*count` is the only output. Otherwise `capacity` must be at least that
 * length or H5C_ERR_SHAPE_MISMATCH is returned; `count` may then be NULL.
 *
 * Prefer h5c_poffset() when you only need your own block. Reach for this when
 * you need every rank's boundaries, for instance to describe or re-create the
 * decomposition.
 *
 * Collective: every rank must call it.
 */
h5c_status_t h5c_ppartition(h5c_file_t *file, const char *path,
                            int64_t *bounds, size_t capacity, size_t *count);

/*
 * Interleaved multi-component write, distributed along the split axis.
 * Each rank supplies `ncomp` component pointers of `n` local elements; the
 * result is an [total_n, ncomp] dataset under path/data. See the pack-versus-
 * stride note in h5c.h: writes pack, reads stride.
 *
 * CALLER BEWARE, and the same applies to the serial calls in h5c.h:
 *
 *   - `n` CANNOT BE VERIFIED. Each component is a bare pointer with no length,
 *     so an `n` larger than a component's real length reads or writes past its
 *     end with no diagnostic. C offers nothing better here; a C++ or Fortran
 *     wrapper should carry lengths with the buffers and check them before
 *     calling. h5cpp passes spans for exactly this reason.
 *   - The zero rules are ASYMMETRIC on purpose: `ncomp == 0` is an error,
 *     because there is no component to write, while `n == 0` is legal and
 *     means this rank owns no rows. With `n == 0` the component pointers
 *     themselves may be NULL.
 *   - `const void *const *comps` forces a cast from `const T *const *` in any
 *     typed wrapper. That is deliberate: a void-pointer array is the only way
 *     C can take a list of differently-typed buffers, and the alternative
 *     would be a per-type entry point for every rank of every type.
 */
h5c_status_t h5c_pwrite_interleaved(h5c_file_t *file, const char *path,
                                    const void *const *comps, size_t ncomp,
                                    size_t n, h5c_type_t type, unsigned flags);

h5c_status_t h5c_pread_interleaved(h5c_file_t *file, const char *path,
                                   void *const *comps, size_t ncomp,
                                   size_t n, h5c_type_t type);

#ifdef __cplusplus
}
#endif

#endif /* H5C_MPI_H */
