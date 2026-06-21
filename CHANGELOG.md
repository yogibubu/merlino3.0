# Changelog

## 2026-06-21 - SEfit/GICForge freeze candidate

This freeze closes the current Merlino semiexperimental-geometry development
cycle relative to the historical MSR workflow.

### Added

- Automatic SEfit diagnostic warning output in `semiexp_warnings.csv`, the text
  report, the HTML report, and the run manifest.
- Warning checks for rank-deficient fits, small weighted-Jacobian singular
  values, ill-conditioned planar component pairs, low robust isotopologue
  weights, and large propagated primitive-coordinate uncertainties.
- Restartable SEfit checkpoints, grouped robust losses, exact
  leave-one-isotopologue-out diagnostics, SVD coordinate-combination diagnostics,
  and constraint audit tables.
- Native semiexperimental input contracts for Cartesian geometries,
  Gaussian-style primitive constraints, isotopologue corrections, QM predicates,
  and parameter classes.
- Legacy MSR compatibility for `.msr` and `.msr.inp` inputs, including frozen
  Z-matrix variables converted to primitive constraints.
- Paper benchmark inputs, generated summaries, and regression snapshot support
  for glycolaldehyde, cyclopentadiene, nitrobenzene, azulene, and norcamphor.

### Changed

- SEfit now treats GIC definitions as frozen schemas during ordinary
  least-squares iterations and rebuilds B matrices only when required.
- Planar molecules use only two independent rotational components, selected by
  rank and stability, with the third component retained as a diagnostic.
- The Fortran77 semiexperimental numerical kernel mirrors the Python scaling,
  rank-revealing linear algebra, and grouped robust-weight logic.
- The GUI exposes the same semiexperimental controls as the command-line
  workflow, including robust loss, restart/checkpoint, leave-one-out, coordinate
  model, constraints, predicates, and parameter classes.

### Validation

- `python3 -m pytest gui/tests -q`
- `git diff --check`
- `latexmk -pdf -interaction=nonstopmode main.tex`
- `latexmk -pdf -interaction=nonstopmode supporting_information.tex`
