# Merlino4 Refactor Plan

## Baseline

Merlino4 starts from the final Merlino3 working baseline:

- Code baseline: `c704017 feat(dvr): add fortran77 path solver`
- Closure docs: `836a418 docs: clarify merlino3 closure point`
- Final Merlino3 freeze: `archives/freezes/merlino3.0_freeze_20260620_141839.tar.gz`

Merlino3 is frozen. Do not continue structural work there.

## Refactor Goal

Merlino4 must make the main functional areas independently maintainable. The GUI
should orchestrate workflows, not contain scientific or backend-specific logic.
Fortran programs should be called through narrow Python wrappers. Gaussian and
DVR parsing/writing should have one source of truth.

## Target Packages

- `merlino_core`: configuration, paths, logging, manifests, job status and common
  exceptions.
- `merlino_geometry`: atoms, XYZ, masses, isotopes, topology, rings and symmetry
  data structures.
- `merlino_gic`: ring numbering, GIC construction, Gaussian GIC coordinate
  definitions and Python/Fortran comparison helpers.
- `merlino_gaussian`: Gaussian input writing, log parsing, scan extraction,
  GIC-value extraction and job metadata.
- `merlino_fortran`: executable discovery, build checks, subprocess wrappers and
  normalized error reporting for GICForge and DVR.
- `merlino_dvr`: Gaussian-log/grid to DVR workflows, Cremer-Pople mapping,
  Fortran bridge integration and output readers.
- `merlino_gui`: PySide6 windows and controllers only; all scientific work goes
  through service interfaces.
- `merlino_data`: local libraries, catalog indexes and curated chemical data.

The current Merlino3 folders stay in place until each area has been migrated
behind tests.

## Interface Contracts

### GICForge

Input:

- Cartesian XYZ file.
- Optional backend options: symmetry, B matrix, Gaussian route/resources.

Output:

- Gaussian `.gjf`.
- Readable GIC report.
- Optional B matrix.
- Manifest JSON with executable, input checksum, options and output paths.

### DVR

Input:

- Completed Gaussian log or prepared grid CSV.
- Boundary/solver settings.
- Optional ring/Cremer-Pople labeling settings.

Output:

- Levels CSV.
- Vectors/profile CSV.
- Summary text.
- Figures where requested.
- Manifest JSON with Gaussian-log checksum, solver, backend and output paths.

### Gaussian

Input:

- Molecular model plus route/resources/workflow options.

Output:

- `.gjf` input.
- Completed `.log` where run locally.
- Parsed scan/path records.
- Manifest JSON.

## Migration Sequence

1. Add the new package skeletons and compatibility imports.
2. Move path/config/logging helpers into `merlino_core`.
3. Consolidate geometry/topology/ring primitives in `merlino_geometry`.
4. Move GIC generation and ring numbering into `merlino_gic`.
5. Introduce `merlino_fortran` wrappers for `gicforge.x` and `path_dvr.x`.
6. Move Gaussian parsing/writing into `merlino_gaussian`.
7. Move DVR workflow orchestration into `merlino_dvr`.
8. Rewire GUI controllers to call service interfaces.
9. Remove compatibility wrappers only after tests cover the new imports.

Each step should end with a small commit and a green validation run.

## Rules

- Do not add new scientific logic directly inside GUI classes.
- Do not duplicate Gaussian log parsers.
- Do not call Fortran executables directly from GUI code.
- Do not hardcode absolute paths except in local environment setup.
- Keep all file-based backend contracts documented and tested.
- Prefer integration tests that exercise real input/output files.
- Keep GICForge and DVR Fortran77 backends independent compiled tools.

## First Milestone

The first real milestone is a no-behavior-change extraction of `merlino_core`
and `merlino_fortran`, with the existing GUI and command-line workflows still
passing `./freeze_check.sh`.
