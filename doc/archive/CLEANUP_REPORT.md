# Cleanup Report (Merlino 4.0)

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
- Removed historical Merlino3-only material from the Merlino4 source tree.
- Updated active docs/UI labels from `Merlino 2.x` to `Merlino 3.0`.

## Residual Duplicates (intentional or low priority)
- Empty `__init__.py` files across packages (benign).
- `geometry/elements.py` and `merlino_fit/topology/elements.py` (currently identical, domain-level split).
- `merlino_fit/topology/readme_topology` and `topology/readme_topology` (compatibility pointers to `doc/topology_doc.txt`).
- Binary/data duplicates:
  - `bin/gicforge.x` and compatibility alias `bin/prova.x`
  - `geometry/fchkin` (sample FCHK kept as canonical; root duplicate removed)

## Historical Material Policy
Merlino4 does not duplicate historical source trees, legacy benchmark scripts,
or old root notes. Use the frozen Merlino3.0 tree for those archives.
