# VPT2/VCI Core

Merlino4 separates Gaussian parsing from the numerical solvers.

## Inputs

- Gaussian FCHK: masses, Cartesian Hessian and Gaussian anharmonic arrays.
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

- `merlino_vpt2_vci.gaussian_qff`: FCHK and normalized-QFF readers.
- `merlino_vpt2_vci.harmonic`: independent Wilson-GF linear algebra.
- `merlino_vpt2_vci.vci`: product-basis VCI matrix elements and dense
  diagonalization for small spaces.
- `merlino_vpt2_vci.workflow`: Gaussian-FCHK to GF/VCI orchestration.

## Fortran77 Backend

- `fortran/vpt2_vci/gf_core.f`: independent GF helper.
- `fortran/vpt2_vci/vci_core.f`: product-basis and dense VCI helpers.

The Fortran code is fixed-form Fortran77 and receives arrays only. It is not a
GDV wrapper.

## Next Numerical Step

The dense VCI path is for small validation spaces. Large VCI spaces need a
standalone Davidson matrix-vector driver using the same VCI matrix-element
contract, with `utilnz.F:NHDiag` used only as a reference source.
