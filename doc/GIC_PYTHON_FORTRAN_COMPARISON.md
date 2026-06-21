# Python/Fortran GIC Identity

Date: 2026-06-20

Updated: 2026-06-21

The production contract is now identity, not compatibility.  GICForge Fortran77
is the canonical Cartesian-to-GIC definition engine.  The Python GIC command
`python -m survibfit.cli gic` runs GICForge and writes the exact ReadGIC lines
from the Fortran `gauin` file.  The legacy pure-Python local builder is retained
only behind `--python-local` for diagnostics and low-level numerical tests.

## Cases Checked

| Case | Python result | Fortran result | Notes |
| --- | --- | --- | --- |
| C4 chain | exact Fortran ReadGIC lines | exact Fortran ReadGIC lines | Regression test compares Python CLI output with `gauin`. |
| OOO triangle | exact Fortran ReadGIC lines | exact Fortran ReadGIC lines | Regression test enforces 3 non-redundant coordinates. |
| Nitrobenzene MSR geometry | exact Fortran ReadGIC lines | exact Fortran ReadGIC lines | Z-matrix dummy conversion now preserves the closed C6 ring; raw and symmetrized GIC counts are 36. |

## Fortran Findings

- Initial GIC generation is the single source of truth: primitives are generated
  by coordinate family, ring coordinates are added before reduction, and final
  redundancy removal is type-local.
- `SYMMALL`/post-check symmetry metadata is part of the frozen GIC definition
  and is propagated to downstream Python programs through `GICDefinition`.
- The Python CLI no longer writes an independently generated ReadGIC basis by
  default, so ring and out-of-plane subspaces cannot diverge from Fortran.

## Current Policy

Use Python for orchestration, schema freezing, B-matrix evaluation, fitting and
post-processing.  Use GICForge for every production Cartesian-to-GIC definition.
The pure-Python local builder is non-authoritative and must not be used for
Python/Fortran comparison or production Gaussian input generation.
