# Semiexperimental Geometry Fortran77 Kernel

`semiexp_core.f` contains independent Fortran77 numerical kernels used to keep
the semiexperimental geometry workflow reproducible outside the Python solver.

Provided routines:

- `M4SEBondB`: analytic Wilson B row for bond distances.
- `M4SEAngleB`: analytic Wilson B row for valence angles.
- `M4SERotConst`: principal moments and rotational constants from Cartesian
  coordinates and isotope masses.
- `M4SENormalEq`: weighted least-squares normal equations, covariance and
  Gauss-Newton Hessian.
- `M4SEClassNormalEq`: the same normal-equation kernel after compressing
  shared parameter classes and blocking fixed classes.

Compile check:

```bash
./compile_check
```
