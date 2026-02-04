# Cleanup Report (Merlino 3.0)

## Completed
- Unified topology source of truth to `merlino_fit/topology`.
- Converted `topology/` into compatibility wrappers.
- Consolidated imports to `merlino_fit.topology` in core/test GUI paths.
- Removed duplicate root helper scripts:
  - `test_provin_writer.py`
  - `test_prova.py`
- Moved root smoke scripts into `scripts/`:
  - `scripts/smoke_geometry_pipeline.py`
  - `scripts/smoke_advanced_window.py`
- Updated references:
  - `readme_rotvib`
  - `FREEZE_NOTES.md`
- Updated active docs/UI labels from `Merlino 2.x` to `Merlino 3.0`.

## Residual Duplicates (intentional or low priority)
- Empty `__init__.py` files across packages (benign).
- `geometry/elements.py` and `merlino_fit/topology/elements.py` (currently identical, domain-level split).
- `merlino_fit/topology/readme_topology` and `topology/readme_topology` (compatibility pointers to `doc/topology_doc.txt`).
- Binary/data duplicates:
  - `bin/prova.x` and `fortran/prova/prova`
  - `fortran/read_xyzin.f` and `fortran/symmetry/read_xyzin.f`
  - `fortran/symmetry/xyzin` and `fortran/symmetry/xyzin2`
  - `h2o.fchk` and `geometry/fchkin`

## Legacy Historical Docs Kept As-Is
- `doc/merlino2.1_freeze.md`
- `doc/addendum_merlino2.1.1`

These intentionally keep original version naming for archive traceability.
