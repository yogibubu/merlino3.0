# Worktree Triage — Merlino 3.0

Date: 2026-03-24

## Goal

Separate active development lines from runtime noise and repository-structure
issues before attempting a deeper cleanup.

## Active Development Lines Detected

### 1. GUI / Gaussian / DOS-Q(T)

Main files:

- `gui/main_window.py`
- `gui/input_panel.py`
- `gui/gaussian.py`
- `gui/project_manager.py`
- `gui/status_reporter.py`
- `gui/xyzin_utils.py`
- `gui/bdpcs3_workflow.py`

Meaning:

- this is the current user-facing Merlino 3.0 workflow line
- it includes Gaussian-property separation, status/reporting, and DOS/Q(T)

### 2. Fragment Pipeline / Delta Correction / Rotational Bridge

Main files:

- `gui/fragment_pipeline_window.py`
- `gui/deltavib_alpha_dialog.py`
- `merlino_fit/survibfit/fragment_pipeline.py`
- `merlino_fit/survibfit/fragment_delta_correction.py`
- `merlino_fit/tests/test_fragment_pipeline.py`
- `merlino_fit/tests/test_fragment_delta_correction.py`

Meaning:

- this is an active chemistry workflow line, not cleanup noise
- the new untracked files should be treated as candidate source files to keep
- the `deltavib_alpha_dialog.py` piece should not be interpreted as the main
  scientific home of the rovibrational work
- that scientific line now lives primarily in `CeDiTT + alpha_resonances`
- in `merlino3.0` it should be treated as a bridge/integration layer

### 3. BDPCS3 / Documentation / Symmetry Work

Main files:

- `merlino_fit/survibfit/modify_geom.py`
- `merlino_fit/survibfit/symmetry_classifier.py`
- `merlino_fit/tests/test_bdpcs3.py`
- `merlino_fit/docs/newbdpcs3.tex`
- `merlino_fit/docs/synthons_bdpcs3_*.tex`
- `merlino_fit/README.md`
- `merlino_fit/pyproject.toml`

Meaning:

- this is active methodological work and manuscript support

## Repository-Structure Issue

### `working/xyzin` is tracked

This is inconsistent with the intended role of `working/` as an ephemeral
runtime workspace.

Current implication:

- `git status` continues to show runtime drift inside a directory that should
  normally stay disposable

Recommendation:

- do not change this automatically during active development
- handle it later as a dedicated git cleanup step:
  - verify no code path relies on a tracked bootstrap file
  - then stop tracking `working/xyzin` explicitly

## Current Safe Cleanup Boundary

Safe now:

- ignore local data-library directories
- document active code lines
- reduce ambiguity around repo structure

Not safe now without review:

- mass-moving files
- pruning modified source files
- changing tracked runtime files such as `working/xyzin`
- collapsing GUI and `merlino_fit` trees without a dedicated integration pass
- expanding the local `DeltaVib/alpha` bridge into a parallel scientific line
