# 2026-09-06: Cleanup decisions and evidence behind the split

This append-only record distinguishes repository evidence from the owner's
account of abandoned attempts. This documentation session ran no builds,
tests, or distributed jobs. Results below are historical commit reports,
not new measurements.

## Public headers and the rejected handle workaround

Commit 5eb61df explains why the previous H5C_HAVE_PARALLEL-dependent public
headers forced separate serial and parallel install prefixes: choosing a
build changed the API exposed by the installed headers. Separating the
parallel declarations made both APIs available from one installation without
replacing its headers. The commit reports serial quick 5/5, parallel quick
5/5, two parallel checks at four ranks, and a serial consumer built with a
normal compiler against the parallel installation's serial target.

An intermediate attempt self-declared hid_t as int64_t so that h5c.h need not
include <hdf5.h> at all, the aim being to keep MPI headers out of the serial
API. It was reverted before any commit, which is why a git history search for
`typedef int64_t hid_t` finds nothing.

Two things settled it, both measured rather than argued:

- hid_t's underlying type is not fixed. This machine carries both forms:
  /opt/system/app/intelpython/2024.2.0/include/H5Ipublic.h has
  `typedef int64_t hid_t;` and
  /opt/system/app/grads/2.0.2/supplibs_rhel8/include/H5Ipublic.h has
  `typedef int hid_t;`. A self-declaration silently disagrees with the second.
  The specific release that changed it was not established here and should not
  be quoted without checking.
- The workaround was unnecessary anyway. A program including <hdf5.h> and
  linking only h5c::h5c_serial builds and runs with a plain compiler against a
  PARALLEL installation, because the imported HDF5 target carries the MPI
  include paths. This was compiled and run, not inferred.

Anyone tempted to repeat the optimization should read those two points first.

## Keeping one visualization body

5eb61df chose five parallel hooks instead of a second writer body. Inspection
of src/serial/h5c_viz.c and src/h5c_viz_internal.h confirms the serial file
holds the writer and a NULL ops pointer selects serial defaults. The five
seams were status agreement (return the local status), transfer properties
(H5P_DEFAULT), empty selections (no parallel select_none hook; serial
zero-row transfers are skipped), counts (zero offsets and local totals), and
tile agreement (leave the local count unchanged). These are fallbacks for an
absent table, not support for missing entries in a non-NULL table. This kept
layout decisions in one implementation while permitting both opening modes.

578676c replaced gather_counts' Allgather, allocated array, and prefix loop
with Exscan for offsets and Allreduce for totals. That removed the allocation
and associated failure/cleanup branches. Rank zero's undefined Exscan result
is explicitly reset. A single rank cannot exercise the inter-rank prefix;
meaningful verification needs more than one rank. The commit reports
passing test_parallel and test_pviz at four ranks, not a timing benchmark.

## A passing suite that missed interoperability

8e51630 documents phase E's absolute path into one developer's home directory.
Other checkouts skipped the real Fortran artifact while the suite still said
"all checks passed"; the named directory had disappeared locally too. Thus
those passes did not exercise actual cross-library interoperability, despite
the project's purpose. The commit reports checking both the missing-artifact
case and a successful run against an artifact from h5fortran's serial test.
Rejecting an unconditional artifact requirement preserved optional sibling
checkouts; explicit skip reporting was chosen to stop hiding the coverage gap.

## Rationale recovered during changelog cleanup

We reduced CHANGELOG from 145 to 32 lines. Git placed the target split at
5eb61df and the serial visualization change at e21ed26. The latter had already
listed pread_rows, but 5eb61df records its implementation. We used the owner's
explicit release boundary, retaining the serial visualization change under
v0.1.0 and the target split and row-read entry under v0.1.1. Banner cleanup,
Exscan, and path repair were not existing changelog bullets; adding them was
rejected under the instruction to reorganize rather than invent entries.

The old changelog explained these choices, which should not be lost:

- Explicit communicator selection and transfer mode avoided hidden parallel
  policy; global validation agreement prevented one rank returning while
  peers entered a collective operation. Empty parallel selections used
  H5Sselect_none to avoid implementation-sensitive zero-length hyperslabs.
- An int8-based bool enum saved space and remained readable by inspection
  tools. The original interoperability claim incorrectly generalized a
  measured enum-to-integer conversion to the reverse direction. The old
  FORMAT text acknowledges this mistake; integer-to-enum reads failed, and
  the read path used native int8 instead. We removed that historical paragraph
  from FORMAT, retaining the current conversion contract.
- Open mode and dataset replacement were separated rather than overloading
  one policy. Enum-based type dispatch avoided repeating the writer per type;
  the old changelog compared ten functions with roughly sixty generated
  Fortran procedures, a historical comparison not recounted here as a fresh
  measurement. Real128 was excluded because the C representation required a
  compiler extension. One-dimensional partition boundaries could not describe
  simultaneous decomposition along multiple axes.
- Product-major coupling was rejected for the shared scheme because an
  unrelated product release must not break the common format. The old
  comparison with h5fortran was historical; its v1.1.1 changelog now records
  independent versions. FORMAT now states only the independent contract.
- 9704297 rejected single-element datasets, repeated scalar attributes, and
  packed strings for short metadata vectors because their call sites were
  awkward. The old changelog recommended datasets for large payloads because
  attributes occupy object headers. The commit reports Fortran reading C
  attributes and C++ reading Fortran attributes, not merely matching schemas.
- Reference HDF5 files were generated rather than committed so expected
  layouts remained inspectable in source. The original four-rank test report
  was a coverage limit, not proof for arbitrary process counts.

The initial cleanup moved rationale into FORMAT/USAGE under the original
instruction. The owner's correction rejected that destination for history.
We moved the attribute-storage rationale and empty-selection explanation
here; scheduler failure reasoning is recorded once in the h5fortran log.
We retained compact current attribute APIs, storage facts, unsupported-feature
limits, and the requirement to create output directories before submission.
The stale "array attributes unsupported" statement was not restored.

## Operational traps supplied by the owner

These two incidents were supplied in the course correction, not reproduced
by this documentation session:

- An argv-based pgrep search can match the shell running that search because
  the shell's command text contains the pattern. Bare pgrep normally searches
  process names; the argv trap concerns full-command-line matching. File
  modification times or scheduler job status should be the primary evidence.
  If such a process search is unavoidable, a pattern like `[w]riter` avoids
  matching its own literal pattern in the search command's argv. A reported
  match alone was therefore rejected as proof of ongoing work.
- The owner reports that this environment's permission filter rejects command
  text containing certain distributed-launcher names, even in instructions
  prohibiting their use. Such an instruction can be refused before delivery.
  Describe the restriction as "no distributed job execution" without spelling
  those names. No filter rejection or bypass was attempted in this session.
