# Final Cleanup Status (Merlino 3.0)

## Final Stop Point
- Final Merlino3.0 commit before starting Merlino4.0: `c704017`
  (`feat(dvr): add fortran77 path solver`).
- Final freeze archive:
  `archives/freezes/merlino3.0_freeze_20260620_141658.tar.gz`.
- Final freeze checksum:
  `81222a3624298d0bd77e51a8b26a6010a5b754471dfd6591058e1d4766687a67`.
- Final validation command: `./freeze_check.sh`.
- Final validation result:
  - GUI tests: 12 passed
  - `merlino_fit` tests: 82 passed
  - `puckering_dvr` tests: 1 passed

Merlino3.0 is closed as a working baseline. New structural work should start in
Merlino4.0 with a full refactor aimed at isolating the main functional areas so
they can evolve independently.

## Merlino4.0 Starting Direction
- Split the current monolithic runtime into independent packages/modules for
  GUI, geometry/topology, GIC/Gaussian generation, Fortran backends, DVR and
  data libraries.
- Define narrow interfaces between Python orchestration and Fortran executables:
  file formats, command-line contracts, output manifests and error reporting.
- Keep GICForge and the Fortran77 DVR as compiled backend tools; keep RDKit,
  SMILES, Gaussian-log parsing, path construction and high-level workflows in
  Python.
- Add integration tests around the interfaces instead of testing only internal
  implementation details.
- Treat `working/`, `projects/` and `archives/freezes/` as local runtime/data
  areas, not source layout examples for Merlino4.0.

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
- References updated (`doc/legacy/readme_rotvib`,
  `doc/legacy/FREEZE_NOTES_2026-02-04.md`).

## Test Stability
- Sign/format-robust regression checks added for GIC output in:
  - `merlino_fit/tests/test_regression_gic.py`
- GUI test dependency installed in venv:
  - `pytest-qt`

## Path DVR / Puckering Integration
- Runtime copy imported under `puckering_dvr/` from
  `/Users/vincenzobarone/Desktop/puckering_dvr_github`.
- The imported copy excludes source `.git`, legacy zip package, and
  bibliography PDFs.
- The main toolbar exposes a dedicated `DVR` window.
- The window runs `puckering_dvr/scripts/mw_path_dvr.py` on any Gaussian
  optimized scan/path log and writes DVR outputs/figures under the selected
  working directories.
- The window can select the newest Gaussian log, run a lightweight preflight,
  launch Gaussian from `gauin.gjf`, chain Gaussian -> DVR, preview DVR results,
  and write a diagnostic `*_dvr_run_manifest.json`.
- Ring puckering labels are optional: Gaussian `QPck/PhiP` values are preserved
  as `gic_*` columns and can be mapped to generalized Cremer-Pople components.
- The intended sequence is:
  `merlino-set` -> `merlino-run` -> generate Gaussian -> run Gaussian -> DVR.
- The Fortran77 DVR backend lives under `fortran/dvr/` and builds
  `bin/path_dvr.x`.
- The GUI exposes Python DVR solvers plus `fortran-sinc-dvr` and
  `fortran-gaussian`.
- The Fortran DVR diagonalizer is `DVRHQRII`, a renamed local copy of HQRII1;
  Jacobi diagonalization is not used.

## Validation Snapshot
- Final `freeze_check.sh` snapshot on 2026-06-20:
  - `gui/tests`: 12 passed
  - `merlino_fit/tests`: 82 passed
  - `puckering_dvr/tests`: 1 passed

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
