# Fortran Backends

Merlino4 keeps only active Fortran backends under `fortran/`. Historical
standalone programs remain in the frozen Merlino3.0 tree.

## Active GICForge Backend

- Source: `fortran/gicforge/`
- Runtime executable used by launchers: `bin/gicforge.x`
- Compatibility runtime alias: `bin/prova.x`
- Build command:

```bash
cd fortran/gicforge
./compile_MAC
```

`compile_MAC` uses `gfortran -std=legacy -w` because the source is intentionally
fixed-form legacy Fortran. This suppresses compiler deprecation noise such as
shared DO labels, arithmetic IF, and Hollerith constants while preserving
failure on actual compile/link errors.

The source-tree build output is `fortran/gicforge/build/gicforge`; runtime
copies are written to `bin/gicforge.x` and `bin/prova.x`. Build logs are written
to `fortran/gicforge/build/gicforge_build.log`. The whole `build/` directory is
ignored by git.

GICForge receives only Cartesian XYZ input from Merlino, builds redundant and
non-redundant GICs, emits a readable report, and creates Gaussian-readable GIC
input. Merlino wraps this as the `gic-define` utility, which writes a frozen
`merlino.gic.definition.v1` schema containing primitives, GIC coefficients,
labels, irreducible representations and the Gaussian block.

B-matrix construction is a separate library contract. The `gic-bmatrix` utility
and `merlino_gic.evaluate_gic_definition` read the frozen schema and evaluate
GIC values and analytic Wilson B rows on the current Cartesian geometry without
rerunning topology perception, redundancy removal or symmetry assignment. GUI
orchestration, RDKit/SMILES, any conversion to Cartesian coordinates, DVR and
post-processing remain Python responsibilities.

The active build intentionally does not compile old Z-matrix, FITPOT/VCI/DVR,
MSR/isotope or rate utilities. Those historical programs are not duplicated in
Merlino4.

## Active DVR Backend

- Source: `fortran/dvr/`
- Runtime executable used by launchers: `bin/path_dvr.x`
- Build command:

```bash
cd fortran/dvr
./compile_MAC
```

The Fortran77 DVR backend reads only `dvrin`. Python remains responsible for
reading Gaussian output, constructing path/grid CSV files and launching the
bridge script. The backend currently supports one-dimensional grid DVR,
one-dimensional distributed Gaussian basis DVR and a two-dimensional product
grid DVR.

Diagonalization is performed by `DVRHQRII` in `fortran/dvr/dvr_hqrii.f`, a
renamed local copy of GICForge `HQRII1`. This avoids Jacobi diagonalization for
large Hamiltonians.

## Merlino4 VPT2/VCI Core

`fortran/vpt2_vci/` now contains new Fortran77 source kernels written for
Merlino4:

- `gf_core.f`: Wilson-GF helper for already independent coordinates.
- `vci_core.f`: small dense VCI helpers and product-basis generation.
- `davidson_core.f`: Davidson support routines independent from Gaussian/GDV.

These routines do not parse external electronic-structure files and do not
build coordinates. Python owns adapter parsing, canonical Merlino input
normalization and workflow orchestration. The Fortran kernels receive numerical
arrays only; Davidson is
implemented with a Merlino4 `matvec + diagonal` contract and is not copied from
Gaussian/GDV.

No historical Gaussian/GDV harmonic Fortran source is kept in Merlino4.

## Semiexperimental Geometry Core

`fortran/semiexp/` contains the Fortran77 numerical kernel for
semiexperimental equilibrium geometry fitting:

- `semiexp_core.f`: analytic B rows for basic internal coordinates, rotational
  constants from Cartesian coordinates and isotope masses, weighted normal
  equations, covariance, least-squares Hessian and the single-substitution
  Kraitchman coordinate kernel using the same reduced substitution mass as the
  Python semiexperimental solver.
- `compile_check`: fixed-form compile check producing only local build
  artifacts.

Python remains responsible for CSV/XYZ parsing, GIC construction, isotope
bookkeeping, line search and workflow manifests. The Fortran source receives
only numerical arrays and is registered as a source backend named `semiexp`.

The `fortran/` root is intentionally documentation-only. Generated compiler
logs such as `error` files are not source and should not be committed.
