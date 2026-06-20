# Repository Layout

This file is the operational map for the Merlino repository.

## Root Policy

The repository root contains only launchers, freeze scripts, git metadata, and
the top-level README. New implementation files should not be added at root.

## Runtime Code

- `gui/`: PySide6 user interface.
- `advanced/`: advanced workflow panel, dedicated DVR window, launchers, and
  viewers.
- `geometry/`: core molecular geometry, rotational/vibrational/thermodynamic
  routines, and `xyzin` utilities.
- `merlino_fit/`: topology, synthons, fitting, BDPCS3 correction, Gaussian GIC
  generation, and test suite.
- `puckering_dvr/`: DVR analysis backend for completed Gaussian outputs.
- `fortran/`: active Fortran backends. Code lives in `fortran/gicforge` and
  `fortran/dvr`; historical standalone utilities remain in Merlino3.0.

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
- Historical Merlino3 notes are not duplicated in Merlino4; use the frozen
  Merlino3.0 tree for that archive.

## Fortran Build Policy

Fortran sources are kept in fixed-form Fortran-compatible style. Build scripts
should use `-std=legacy` and suppress compiler deprecation noise from historical
constructs, while still failing on real compilation/link errors.

`fortran/gicforge/compile_MAC` is the canonical GICForge build command. It updates:

- `fortran/gicforge/build/gicforge`
- `bin/gicforge.x`
- `bin/prova.x`

`fortran/dvr/compile_MAC` is the canonical Fortran DVR build command. It updates:

- `fortran/dvr/build/path_dvr`
- `bin/path_dvr.x`
