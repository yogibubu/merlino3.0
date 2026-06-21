# GF And VPT2/VCI Core

Merlino4 separates external file parsing from the numerical solvers.

## Inputs

- `merlino_gf.HessianInput`: canonical Merlino Cartesian geometry, masses and
  Hessian for the harmonic GF branch.
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

- `merlino_gf.models`: canonical Cartesian Hessian input model for GF/PED.
- `merlino_gf.harmonic`: independent Wilson-GF linear algebra.
- `merlino_gf.internal`: Cartesian Hessian plus Merlino frozen GIC/B matrix to
  GF frequencies, normal modes and PED. It supports both on-the-fly Merlino GIC
  construction and frozen `merlino.gic.definition.v1` schemas from
  `gic-define`.
- `merlino_gf.service`: GF/PED report and CSV service. The CSV tables retain
  GIC names and irreps from the frozen schema.
- `merlino_vpt2_vci.models`: canonical Merlino anharmonic input data models.
- `merlino_vpt2_vci.gaussian_qff`: Gaussian FCHK adapter and normalized-QFF
  reader.
- `merlino_vpt2_vci.vci`: product-basis VCI matrix elements and dense
  diagonalization for small spaces.
- `merlino_vpt2_vci.vpt2`: VPT2 energies and VPT2/VCI comparison on the same
  canonical QFF and mode-selection options.
- `merlino_vpt2_vci.davidson`: independent symmetric Davidson diagonalizer
  using only `matvec` and an approximate diagonal.
- `merlino_vpt2_vci.workflow`: compatibility Gaussian-FCHK to VCI orchestration.

## Fortran77 Backend

- `fortran/vpt2_vci/gf_core.f`: independent GF helper.
- `fortran/vpt2_vci/vci_core.f`: product-basis and dense VCI helpers.
- `fortran/vpt2_vci/vpt2_core.f`: VPT2 helper for quartic first-order and cubic
  second-order corrections on a supplied basis.
- `fortran/vpt2_vci/davidson_core.f`: independent Davidson support routines.

The Fortran code is fixed-form Fortran77 and receives arrays only. It is not a
GDV wrapper.

Historical Gaussian/GDV Fortran sources are not part of the active Merlino4
tree.

## GF/PED From Cartesian Hessian

GF/PED is a separate harmonic workflow implemented in `merlino_gf`. It may use
GICs and B matrices; it is not part of the VPT2/VCI package.

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

The production frozen-coordinate path is:

```text
gic-define reference geometry
-> frozen GIC definition JSON
-> Gaussian FCHK adapter, or another future Hessian adapter
-> evaluate the frozen GIC definition and B matrix on the Hessian geometry
-> transform Hcart to F(GIC)
-> optional Pulay scaling: F_ij <- F_ij sqrt(s_i s_j)
-> Wilson GF frequencies, normal modes and PED
```

For this Hessian transformation, the B matrix is evaluated with Cartesian
coordinates in bohr, consistent with the canonical Hessian units Eh/bohr^2.
The standalone `gic-bmatrix` utility remains a geometry/fit diagnostic and
reports derivatives in the coordinate units of the supplied geometry file.

The `gic-gf` CLI exposes this branch:

```bash
python -m merlino gic-gf \
  --schema gic_definition.json \
  --fchk gauin.fchk \
  --scale-file pulay_scale.txt \
  --out gic_gf_ped_report.txt \
  --csv-dir gic_gf_csv
```

GIC symmetrization is optional and is already decided when the frozen
definition reaches this step. GF evaluates the frozen schema, carries GIC names
and irreps through reports/CSV outputs, and never changes coordinate symmetry.

## VPT2/VCI In Cartesian Normal Modes

The anharmonic branch consumes a quartic force field expressed in Cartesian
normal modes. It does not build, symmetrize or evaluate GICs. Mode selection,
frequency windows, pruning and block separation operate on normal-mode indices
and optional mode symmetry labels supplied with the QFF workflow.

## VCI Basis Control

`VCIOptions` controls the numerical basis before diagonalization:

- active-mode selection and frequency windows
- harmonic basis energy cutoff and maximum basis size
- minimum/maximum quanta for each mode
- minimum/maximum total quanta for one-, two-, three- and four-mode excitation
  classes
- cubic/quartic force-constant pruning
- optional block separation by mode symmetry labels
- coefficient threshold for reporting dominant basis-state contributions

The VCI result reports excitation energies, solved block metadata, expectation
values of modal quanta for each final state, and dominant basis-state weights.

## VPT2/VCI Comparison

`compare_vpt2_vci` applies the same reduced-mode selection and pruning options
to VPT2 and VCI, then reports absolute and excitation-energy differences. VPT2
uses quartic terms at first order and cubic terms at second order in the same
dimensionless normal-coordinate convention used by the VCI Hamiltonian.

Large VCI spaces use the standalone Davidson contract. Gaussian/GDV routines
are not dependencies or copy sources for Merlino4.
