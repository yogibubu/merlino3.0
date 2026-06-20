# Final Cleanup Status (Merlino 3.0)

## Environment and Commands
- `merlino-set` points to `~/merlino3.0` and activates the first available
  Merlino environment: `~/.venvs/merlino`, repo-local `.venv`, or conda
  `$MERLINO_CONDA_ENV` (default `merlino_26`, fallback `merlino`).
- `merlino-run` uses active venv and runs `app.py` from `MERLINO_HOME`.
- `merlino-run-bg` runs GUI in background with log output.
- `merlino-run-check` verifies the GUI runtime dependencies, including
  `PySide6`, `PIL`, `rdkit`, `numpy`, `scipy`, and `matplotlib`.
- Canonical shell sequence:
  `source ~/.bashrc && merlino-set && merlino-run-check && merlino-run`.
- If `merlino-run-check` reports missing GUI packages in the active env, use
  `merlino-install-gui-deps`.
- `merlino-test-all` runs full test suite:
  - `merlino_fit/tests`
  - `gui/tests` (offscreen Qt)

## Structural Cleanup Completed
- Topology source of truth consolidated to `merlino_fit/topology`.
- `topology/` converted to compatibility wrappers.
- Imports in active code/tests migrated to canonical `merlino_fit.topology` paths.
- Root duplicate helper scripts removed:
  - `test_provin_writer.py`
  - `test_prova.py`
- Root smoke scripts moved under `scripts/`:
  - `scripts/smoke_geometry_pipeline.py`
  - `scripts/smoke_advanced_window.py`
- References updated (`readme_rotvib`, `FREEZE_NOTES.md`).

## Test Stability
- Sign/format-robust regression checks added for GIC output in:
  - `merlino_fit/tests/test_regression_gic.py`
- GUI test dependency installed in venv:
  - `pytest-qt`

## Puckering DVR Integration
- Runtime copy imported under `puckering_dvr/` from
  `/Users/vincenzobarone/Desktop/puckering_dvr_github`.
- The imported copy excludes source `.git`, legacy zip package, and
  bibliography PDFs.
- Advanced GUI now exposes `Puckering DVR – Gaussian scan analysis`.
- The panel runs `puckering_dvr/scripts/mw_path_dvr.py` on a Gaussian log and
  writes DVR outputs/figures under the selected working directories.
- The intended sequence is:
  `merlino-set` -> `merlino-run` -> generate/run Gaussian -> run Puckering DVR.

## Validation Snapshot
- `merlino_fit/tests`: 60 passed
- `gui/tests`: 8 passed
- Combined via `merlino-test-all`: all green

## Residual Legacy (Intentional)
- Historical docs kept with 2.1 naming for archive traceability:
  - `doc/merlino2.1_freeze.md`
  - `doc/addendum_merlino2.1.1`
- Low-priority duplicate binaries/data listed in `doc/CLEANUP_REPORT.md`.

## Current Status
Project is ready for feature development on a stable baseline.

## Stable Data Paths
- Semi-experimental structure library (SE): `projects/se_library`
  - copied from: `/Users/vincenzobarone/Downloads/SE`
- PCS2 structure library: `projects/pcs2_library`
  - copied from: `/Users/vincenzobarone/Downloads/PCS2`
