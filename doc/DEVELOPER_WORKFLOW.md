# Merlino4 Developer Workflow

Every new workflow should be added in this order:

1. Define a non-Qt service API with dataclasses for normalized inputs/outputs.
2. Add CLI access under `python -m merlino`.
3. Write a `merlino.run.v1` manifest with input/output checksums.
4. Add focused tests using small fixtures.
5. Add GUI wiring only after the service and CLI are tested.
6. Document file contracts and benchmark commands.
7. Run `./freeze_check.sh` before committing.

Rules:

- GUI classes must not contain scientific algorithms.
- Gaussian, Fortran and other external formats are adapters, not internal data
  models.
- Prefer typed Merlino errors from `merlino_core.errors`.
- Important numerical fits should use the shared `merlino_core.numerics`
  primitives for damped normal equations, step limiting and rank/condition
  diagnostics unless a workflow has a documented reason to use a specialized
  solver.
- Store new project outputs under `inputs/`, `runs/`, `outputs/`, `reports/`,
  `cache/` or `logs/`.
