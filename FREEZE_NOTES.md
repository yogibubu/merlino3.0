# Merlino 3.0 Freeze Notes

- Freeze date: `2026-02-04 19:43:23 CET`
- Workspace root: `/Users/vincenzobarone/merlino3.0`
- VCS status: no `.git` repository detected in this workspace

## Environment

- Python: `Python 3.9.6`
- Pytest: `pytest 9.0.2`

## Verification Commands

Run from workspace root:

```bash
QT_QPA_PLATFORM=offscreen pytest -q gui/tests
PYTHONPATH=/Users/vincenzobarone/merlino3.0/merlino_fit pytest -q merlino_fit/tests
```

Expected at freeze time:

- `gui/tests`: `8 passed`
- `merlino_fit/tests`: `59 passed`

## Notes

- `scripts/smoke_geometry_pipeline.py` and `scripts/smoke_advanced_window.py` are smoke scripts, not pytest test modules.
- Recommended next step after freeze: structural refactor in dedicated branch/snapshot, keeping this freeze as baseline.
