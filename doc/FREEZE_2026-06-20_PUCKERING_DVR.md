# Freeze 2026-06-20 - Puckering DVR Integration

## Scope

This checkpoint integrates the path-DVR backend into the Merlino 3.0 runtime
flow. The DVR backend is generic for Gaussian scans; puckering/Cremer-Pople
handling is an optional labeling and post-processing layer.

Included changes:
- ring-puckering Gaussian GIC generation with inactive `RPck....` and active
  `QPck....`/`PhiP....`
- Fortran `prova` output support for the same `QPck/PhiP` functional GICs
- Python arbitrary-ring GIC generation, tested through seven-membered rings
- vendored runtime copy under `puckering_dvr/`
- Advanced GUI launcher for DVR analysis of Gaussian scan/path logs
- Gaussian GIC value parsing into `gic_*` output columns
- optional fitted bridge from Gaussian `QPck/PhiP` components to generalized
  Cremer-Pople components
- environment checks updated for `scipy` and `matplotlib`

## Environment

Canonical startup sequence:

```bash
source ~/.bashrc
merlino-set
merlino-run-check
merlino-run
```

`merlino-run-check` must report:
- `PySide6`
- `PIL`
- `rdkit`
- `numpy`
- `scipy`
- `matplotlib`

## Puckering Workflow

1. Build or load the molecule in Merlino.
2. Open Advanced Calculations.
3. Generate the Gaussian input with ring GICs.
4. Run Gaussian.
5. Run `Path DVR – Gaussian scan analysis` on the resulting log.

For non-puckering Gaussian scans, leave Cremer-Pople labeling disabled. The DVR
still uses the mass-weighted Cartesian path and does not require any ring
coordinate.

The command-line equivalent is:

```bash
python puckering_dvr/scripts/mw_path_dvr.py \
  --gaussian-log working/gauin.log \
  --log-selection last-per-link \
  --boundary periodic \
  --solver fourier \
  --compute-rotconst \
  --label-cremer-pople \
  --outdir working/puckering_dvr_outputs \
  --figdir working/puckering_dvr_figs \
  --prefix puckering_dvr
```

## Vendored DVR Copy

Source used for the runtime copy:

```text
/Users/vincenzobarone/Desktop/puckering_dvr_github
```

Excluded from the Merlino copy:
- source `.git`
- `puckering_dvr_package_legacy.zip`
- bibliography PDFs

These files are not required for runtime execution or tests.

## Verification

Run from the repository root:

```bash
source ~/.bashrc
merlino-run-check
python -m py_compile advanced/advanced_window.py advanced/launchers/puckering_dvr_launcher.py
python -m py_compile merlino_fit/survibfit/puckering_gaussian.py merlino_fit/survibfit/cli.py
python puckering_dvr/scripts/mw_path_dvr.py --help
PYTHONPATH=merlino_fit python -m pytest -q merlino_fit/tests/test_puckering_gaussian.py
```
