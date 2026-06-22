# Semiexperimental Geometry Fortran77 Kernel

`semiexp_core.f` contains independent Fortran77 numerical kernels used to keep
the semiexperimental geometry workflow reproducible outside the Python solver.
The production workflow is Python-orchestrated because input handling,
automatic GIC generation, symmetry, isotopologue metadata, reporting and GUI
state are managed more robustly there. This Fortran77 layer is intentionally a
validated kernel layer, not a second full workflow implementation.

Provided routines:

- `M4SEBondB`: analytic Wilson B row for bond distances.
- `M4SEAngleB`: analytic Wilson B row for valence angles.
- `M4SERotConst`: principal moments and rotational constants from Cartesian
  coordinates and isotope masses.
- `M4SETrustNormalEq`: weighted SVD-equivalent More-Hebden trust-region
  Levenberg-Marquardt step, covariance and Gauss-Newton Hessian. The SVD
  subproblem is solved in Fortran77 by diagonalizing the scaled symmetric Gram
  matrix, giving the same mathematical step as the Python SVD kernel.
- `M4SENormalEq`: backward-compatible wrapper with inactive trust radius,
  matching the Python rank-revealing LM step without an active radius.
- `M4SEClassTrustNormalEq`: trust-region kernel after compressing shared
  parameter classes and blocking fixed classes.
- `M4SEClassNormalEq`: backward-compatible class wrapper with inactive trust
  radius.
- `M4SERobustGroupWeights`: grouped robust weights matching the Python
  isotopologue-block robust loss convention.

Compile check:

```bash
./compile_check
```
