# User documentation cleanup

Scope: documentation only; no build, tests, jobs, or commits.

- Reduced README to prerequisites, build/install, one C example, and user links.
- Kept API calls and ownership in USAGE; made FORMAT the single reference for
  dimensions, stored types, partition layout, and visualization field ordering.
- Split visualization into USAGE-visualization because mesh setup and XDMF
  postprocessing form a separate user workflow.
- Removed implementation rationale, internal source locations, test procedures,
  and stale claims about h5cpp reading partition data itself.
- Corrected the close-after-use example, C type coverage, and transfer-mode wording.
- Moved docs/TODO.md to log/TODO.md without changing its historical contents.
- Moved the completed visualization brief unchanged to
  log/2026-09-05-visualization-writer.md. The date is the last modification date
  from git log for docs/VISUALIZATION-WRITER.md: e21ed26, 2026-09-05.
- The brief contains no unique on-disk contract missing from FORMAT. Its input
  validation requirements belong to the visualization usage guide. Its old
  source paths and requested work remain historical, not current instructions.
- Checked local Markdown links and heading anchors, whitespace, and changed paths.
  Confirmed both relocated records match their original Git contents byte for byte.

Line counts: README 253 -> 54; FORMAT 263 -> 111; USAGE 489 -> 325;
USAGE-visualization new -> 90; TODO moved 69 -> 69; visualization brief moved 42 -> 42.
