# Finiteness of Planar Central Configurations of Five Bodies

This repository contains the LaTeX manuscript, its compiled PDF, and
independent computational checkers reconstructed from the manuscript for
Zhe Chen's draft paper, "Finiteness of planar central configurations of
five bodies for all positive masses."

## Manuscript

- `main.tex`: paper entry point and shared definitions.
- `sections/d1.tex` through `sections/d8.tex`: the eight diagram sections.
- `main.pdf`: inspected 57-page draft, compiled on 5 October 2026.
- `build.ps1`: Windows build script for an existing pdfLaTeX installation.

Compile `main.tex` with pdfLaTeX three times to resolve references.
The source uses standard AMS LaTeX packages.

## Independent computational verification

Requires Python 3.10 or later. From the repository root, run:

```text
python verification/verify_all.py
```

No package installation or network access is needed. Unmodified SymPy
1.14.0 and mpmath 1.3.0 wheels are bundled in `verification/vendor/`;
their hashes are checked before loading. The run takes several minutes
and regenerates the exact execution records.

The recorded 5 October 2026 run passed all **78 computational checks**:
**42 claim groups PASS, 1 not applicable, 0 unfinished, 0 failed**, exit 0.
The two supplementary programs for Remark 11.10 are included and also
passed a separate run from a clean extracted delivery directory:

- `ak_diagram8_ray_diagrams.py`: all ten distance/four-product orders.
- `ak_diagram8_pattern_boundaries.py`: sharp ray window and boundary patterns.

`verification/program_index.json` maps all 54 Table 5 names: 42 directly
implemented programs and 12 documented replacements. These are reconstructed
checkers, not the original author's code. Some original computations are
covered by stronger exact replacements, as described in the reviewer guide.

## Reviewer records

- [Running instructions and coverage](verification/README.md)
- [Reviewer guide and mathematical scope](verification/REVIEWER_GUIDE.md)
- [Authoritative full-run receipt](verification/certificates/verification_run.json)
- [Release summary](verification/RELEASE_SUMMARY.json)
- [Paper program index](verification/program_index.json)

The receipt records the executed check identifiers and source, code and
certificate hashes. A hash identifies the verified revision; mathematical
correctness comes from replaying the checks and reviewing their models.
Computational coverage does not formalize the analytic proof or the full
main theorem in a theorem prover.
