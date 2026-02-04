# Implementation Notes

## Analytic vs numerical derivatives
- Bond/angle/dihedral/oop/linear gradients are analytic (dihedral via AD).
- Hessian transforms are available in `survibfit/transforms.py`.
- Curvature term is included via finite-difference in internal space when requested.
- G and its derivatives are computed via finite differences on B.

These are intended for prototyping; analytic implementations can replace the
finite-difference parts later.

## Recent symmetry updates
- Global symmetry detection now uses cached candidate operations, symbol grouping,
  radial pre-filters, and early exits for faster matching.
- Optional inertia-based reduction of `max_n` avoids unnecessary high-order axes.
- Relative tolerance term (`symmetry_tol_rel`) improves robustness for large systems.
- Pattern reports include ring clusters, condensed-ring reductions, and dihedral selection.
