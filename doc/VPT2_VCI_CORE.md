# VPT2/VCI Core

Merlino4 separates external file parsing from the numerical solvers.

## Inputs

- `HessianInput`: canonical Merlino Cartesian geometry, masses and Hessian.
- `AnharmonicInput`: canonical Merlino normal-coordinate anharmonic data.
- Gaussian FCHK: supported only through an adapter that populates canonical
  Merlino inputs.
- Normalized QFF text: optional cubic/quartic normal-coordinate force constants.

The normalized QFF format is intentionally simple:

```text
FREQ 1 100.0
CUBIC 1 1 2 -3.5
QUARTIC 1 1 1 1 4.0
```

Indices are one-based in files and are converted to zero-based indices in
Python. Coefficients are currently interpreted in cm^-1 in dimensionless normal
coordinates.

## Python Backend

- `merlino_vpt2_vci.models`: canonical Merlino input data models.
- `merlino_vpt2_vci.gaussian_qff`: Gaussian FCHK adapter and normalized-QFF
  reader.
- `merlino_vpt2_vci.harmonic`: independent Wilson-GF linear algebra.
- `merlino_vpt2_vci.internal_gf`: Cartesian Hessian plus Merlino
  non-redundant GIC/B matrix to GF frequencies and PED.
- `merlino_vpt2_vci.vci`: product-basis VCI matrix elements and dense
  diagonalization for small spaces.
- `merlino_vpt2_vci.davidson`: independent symmetric Davidson diagonalizer
  using only `matvec` and an approximate diagonal.
- `merlino_vpt2_vci.workflow`: Gaussian-FCHK to GF/VCI orchestration.

## Fortran77 Backend

- `fortran/vpt2_vci/gf_core.f`: independent GF helper.
- `fortran/vpt2_vci/vci_core.f`: product-basis and dense VCI helpers.
- `fortran/vpt2_vci/davidson_core.f`: independent Davidson support routines.

The Fortran code is fixed-form Fortran77 and receives arrays only. It is not a
GDV wrapper.

Historical Gaussian/GDV Fortran sources are not part of the active Merlino4
tree.

## GF/PED From Cartesian Hessian

The tested harmonic path is:

```text
Gaussian FCHK adapter, or another future adapter
-> Merlino HessianInput
-> Merlino topology primitives
-> Merlino non-redundant GIC transform U
-> Bq = U^T B
-> G = Bq M^-1 Bq^T
-> F = A^T Hcart A, A = M^-1 Bq^T G^-1
-> Wilson GF frequencies and PED
```

Gaussian is only the current source adapter in this test. The solver-facing
input is `HessianInput`; the GICs, B matrix, GF transformation and PED are
computed by Merlino.

## Next Numerical Step

The dense VCI path is for small validation spaces. Large VCI spaces use the
standalone Davidson contract. Gaussian/GDV routines are not dependencies or
copy sources for Merlino4.
