# Repository Layout

This file is the operational map for Merlino 3.0.

## Root Policy

The repository root contains only launchers, freeze scripts, git metadata, and
the top-level README. New implementation files should not be added at root.

## Runtime Code

- `gui/`: PySide6 user interface.
- `advanced/`: advanced workflow panel, launchers, and viewers.
- `geometry/`: core molecular geometry, rotational/vibrational/thermodynamic
  routines, and `xyzin` utilities.
- `merlino_fit/`: topology, synthons, fitting, BDPCS3 correction, Gaussian GIC
  generation, and test suite.
- `puckering_dvr/`: DVR analysis backend for completed Gaussian outputs.
- `fortran/`: Fortran backends. Active PROVA is `fortran/prova`; other
  subfolders are specialized or legacy numerical modules.

## Generated And Local Files

- `working/`: GUI/runtime work area. Ignored by git.
- `projects/`: local data libraries. Ignored by git.
- `archives/freezes/`: local freeze archives and checksums. Ignored by git.
- `fortran/*/build/`: compiler logs and build intermediates. Ignored by git.

## Documentation

- `doc/PROJECT_STATUS.md`: current project state.
- `doc/FINAL_CLEANUP_STATUS.md`: cleanup and validation history.
- `doc/RING_NUMBERING_CONVENTION.md`: ring numbering convention shared by
  Python and Fortran.
- `doc/reports/`: larger reports and TeX artifacts.
- `doc/legacy/`: historical root notes retained for traceability.

## Fortran Build Policy

Fortran sources are kept in fixed-form Fortran-compatible style. Build scripts
should use `-std=legacy` and suppress compiler deprecation noise from historical
constructs, while still failing on real compilation/link errors.

`fortran/prova/compile_MAC` is the canonical PROVA build command. It updates:

- `fortran/prova/prova`
- `bin/prova.x`

