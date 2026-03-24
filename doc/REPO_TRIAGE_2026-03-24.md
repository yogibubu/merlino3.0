# Repo Triage — Merlino 3.0

Date: 2026-03-24

## Purpose

This note fixes the current interpretation of the repository layout so that
cleanup and future refactors do not mix active code, local workspaces, and
data libraries.

## Main Active Code Areas

- `gui/`
  - main GUI workflow
  - Gaussian/property input separation
  - DOS/Q(T), fragment pipeline, status/report logic
  - transitional `DeltaVib/alpha` bridge for compatibility with external rovibrational work
- `geometry/`
  - rovibrational and thermo pipelines
  - Gaussian parsing and rotational/vibrational processing
- `merlino_fit/survibfit/`
  - fitting, fragment pipeline, delta correction, geometry-edit logic
- `merlino_fit/topology/`
  - canonical topology implementation

## Supporting Areas

- `advanced/`
  - secondary viewers and advanced windows
- `scripts/`
  - smoke tests and helper scripts
- `doc/`
  - architecture, workflow, cleanup, and testing notes

## Runtime / Ephemeral Areas

- `working/`
  - runtime workspace
  - should not be treated as a source directory
- `.pytest_cache/`
- `__pycache__/`

## Data Library Areas

The top-level `projects/` directory in this repository is not a project registry.
It currently hosts molecule libraries used by the fragment pipeline:

- `projects/se_library/`
- `projects/pcs2_library/`
- `projects/hpcs2_library/`

These should be treated as local data/library content, not as source code.

## Current Rule

When reorganizing `merlino3.0`:

- do not interpret `projects/` using the `brain` meaning
- keep code cleanup separate from data-library management
- prefer changes that do not disturb the active GUI / `merlino_fit` work already in progress
- do not treat the local `DeltaVib/alpha` code as the primary scientific line
- treat `CeDiTT + alpha_resonances` as the primary line for rovibrational / `DeltaVib` methodology
- keep only the compatibility layer that Merlino will need downstream

## Immediate Practical Consequence

`projects/hpcs2_library/` is now aligned with the existing ignore policy used
for `projects/se_library/` and `projects/pcs2_library/`, so the worktree does
not get flooded by untracked library files.
