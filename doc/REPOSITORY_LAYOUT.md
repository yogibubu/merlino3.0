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
- `merlino_gic/`: frozen GIC definition and reusable B-matrix evaluation.
- `merlino_gf/`: harmonic Cartesian-Hessian/GIC-B/GF/PED package, physically
  separated from VPT2/VCI.
- `merlino_vpt2_vci/`: normal-mode anharmonic QFF, VPT2, VCI and Davidson.
- `puckering_dvr/`: DVR analysis backend for completed Gaussian outputs.
- `fortran/`: active Fortran backends. Executable code lives in
  `fortran/gicforge` and `fortran/dvr`; independent source kernels for
  harmonic GF, normal-mode VPT2/VCI and Davidson live in `fortran/vpt2_vci`.

## Generated And Local Files

- `working/`: GUI/runtime work area. Ignored by git.
- `projects/`: local data libraries. Ignored by git.
- `archives/freezes/`: local freeze archives and checksums. Ignored by git.
- `fortran/*/build/`: compiler logs and build intermediates. Ignored by git.

## Documentation

- `doc/PROJECT_STATUS.md`: current project state.
- `doc/PACKAGE_ARCHITECTURE.md`: package boundaries and workflow data flow.
- `doc/RING_NUMBERING_CONVENTION.md`: ring numbering convention shared by
  Python and Fortran.
- `doc/archive/`: cleanup, freeze and triage history.
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

`fortran/vpt2_vci/compile_check` verifies the independent GF/VCI/Davidson
source kernels:

- `fortran/vpt2_vci/build/gf_core.o`
- `fortran/vpt2_vci/build/vci_core.o`
- `fortran/vpt2_vci/build/davidson_core.o`
