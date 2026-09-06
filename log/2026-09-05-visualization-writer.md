# h5c visualization writer implementation requirements

## Purpose

The goal was to enable h5c to output the same HDF5 scheme 1 as `h5fortran`'s `t_phdf5_writer`.
XDMF/XML generation was left to the independent tool `h5xdmf`, without implementing it here.

## References

- `h5fortran/src/fypp/parallel/h5fort_parallel_visualization.fypp`
- `h5fortran/docs/USAGE-visualization.md`
- `h5xdmf/src/h5xdmf/schemes/v1.py`
- `h5xdmf/docs/design.md`

The requirement was compatible on-file group, dataset, attributes, shape, and dtype,
not a mechanical translation of the `h5fortran` API. No build dependencies or submodule for h5xdmf or h5fortran were to be added.

## Required Conditions

- `scheme_version` was required to be the constant `1`, independent of h5c product SemVer.
- The root attributes `scheme_version` and `time` were required in the output.
- Each mesh was to output geometry, optional connectivity, and point/cell data.
- `topology_type`, `nodes_per_element`, and `attribute_type` were to be assigned according to scheme 1.
- Vector/Tensor were to reuse the existing `h5c_write_interleaved()` and be stored as `[n, ncomp]`.
- `ncomp` was to be 1, 3, 6, or 9. Tensor6 was to use the ParaView order `XX, YY, ZZ, XY, YZ, XZ`.
- Parallel was to accept rank-local connectivity as 0-origin and add the global node offset.
- Out-of-range connectivity, inconsistent attributes, and metadata mismatches where rank agreement was required were to fail before writing.
- A small C API that h5cpp could wrap was required, without a duplicate implementation for h5cpp.

## Implementation Approach

The approach was to combine the existing file, attribute, parallel, and interleaved API and add only
visualization-specific code. No new library dependencies, plugin mechanism, configuration files, or XDMF generation were to be added.

## Verification

The requirement was to add 1 test generating a small scheme 1 file and inspect group, attributes, shape,
and values with the HDF5 C API. If possible, verification with `h5xdmf validate` on PATH was to be added as an optional test.
MPI tests were required to run through existing job scripts, with no direct execution of `mpiexec`.

After the changes, README, `docs/USAGE.md`, `docs/FORMAT.md`, `docs/TODO.md`, and
`log/CHANGELOG.md` were to be updated only as needed.
