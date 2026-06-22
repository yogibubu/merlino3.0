# Freeze 2026-06-20 - DVR Analysis-Only Backend

This checkpoint removes Gaussian input/path generation from the vendored
`puckering_dvr` backend.

## Architecture

- Merlino/`merlino_fit` is the only supported place where Gaussian puckering
  paths and `.gjf` files are generated.
- `puckering_dvr/scripts/mw_path_dvr.py` consumes completed Gaussian outputs and
  computes DVR levels, wavefunctions, properties, and optional Cremer-Pople/GIC
  labels.
- The DVR backend remains generic for any optimized Gaussian scan/path; the
  scanned coordinate does not have to be puckering.

## Removed From DVR

- CLI options for Gaussian preparation, including `--prepare-gaussian`,
  `--gjf-out`, `--manifest-out`, phase-grid options, start-pucker options, and
  Gaussian route/resource options.
- Internal Gaussian `.gjf`/manifest writer routines from
  `puckering_dvr/scripts/mw_path_dvr.py`.
- Smoke-test expectations that the DVR creates Gaussian input files.

## Validation

Run on 2026-06-20:

```bash
python -m py_compile puckering_dvr/scripts/mw_path_dvr.py
PYTHON=python bash puckering_dvr/scripts/smoke_test.sh
PYTHON=python ./freeze_check.sh
```

`freeze_check.sh` results:

- GUI tests: 8 passed
- `merlino_fit` tests: 79 passed
- `puckering_dvr` tests: 1 passed

