# Merlino Semiexperimental File Formats

This document defines the Merlino4 semiexperimental equilibrium-geometry file
contracts. The recommended production setup uses two input files:

1. A job file, `*.mse.toml`, containing calculation keywords, Cartesian parent
   geometry and optional fixed-parameter definitions.
2. An isotopologue observations file, `*.toml`, `*.json` or `*.csv`, containing
   ground-state rotational constants, vibrational corrections, electronic
   corrections, isotope substitutions and optional experimental uncertainties.

Gaussian-style Cartesian `.com/.gjf` files remain accepted as geometry inputs
for interoperability. Gaussian Z-matrices and MSR-style variable blocks are not
valid Merlino4 SE geometry inputs.

## 1. Job File

Recommended extension: `.mse.toml`.

Schema identifier:

```toml
schema = "merlino.semiexp.job.v1"
```

Complete minimal example:

```toml
schema = "merlino.semiexp.job.v1"
title = "cyclopentadiene SE fit"

[files]
observations = "cyclopentadiene_isotopologues.toml"

[fit]
backend = "python"
coordinate_model = "gic"
observable = "moments"
rotational_components = "auto"
max_step = 0.25
damping = 1.0e-8
step = 1.0e-4
prune_condition = 0.0

[geometry]
units = "angstrom"
atoms = [
  ["C",  0.0000000000,  1.2109940500, -0.0000437769],
  ["C",  1.1696665200,  0.3421049300,  0.0000360217],
  ["C", -1.1696664900,  0.3421050500,  0.0000319263],
]

[constraints]
fix_hydrogen_parameters = true
fixed_gic_patterns = [
  "bond(1,2)",
]
modredundant = [
  "A 2 1 3 F",
]

[[qm_predicates]]
pattern = "GIC001"
value = 1.234
sigma = 0.010
source = "QM reference"

[[parameter_classes]]
name = "CH_stretches"
mode = "shared"
patterns = ["bond(1,6)", "bond(2,7)"]

[[parameter_classes]]
name = "XYH_angles"
mode = "fixed"
patterns = ["angle"]
```

### Job Tables

`[files]`

- `observations`: path to the isotopologue observations file. Relative paths
  are resolved from the job-file directory.

`[fit]`

- `backend`: `python` or `fortran77`. The GUI asks for this choice per run.
- `coordinate_model`: `gic` or `cartesian_symmetry`. `gic` builds
  non-redundant symmetry-adapted GICs and uses their B matrix.
  `cartesian_symmetry` does not use a Hessian or a B matrix: it removes
  Cartesian translations/rotations from the parent geometry, projects the
  remaining displacement space by the detected point-group irreps and optimizes
  only the totally symmetric Cartesian directions.
- `observable`: `moments`, `rotational_constants` or `auto`. The Merlino
  standard is `moments`.
- `rotational_components`: `auto`, `ABC`, `AB`, `AC` or `BC`. This only matters
  for direct rotational-constant fits.
- `max_iter`: optional explicit iteration cap. If omitted, Merlino uses an
  automatic cap proportional to the number of effective fitted parameters.
- `step`: finite working-coordinate step used only by fallback numerical
  derivatives. Current GIC and symmetry-Cartesian SEfit paths use
  analytic Cartesian derivatives for rotational observables.
- `damping`: initial Levenberg-Marquardt damping for the trust-region solver.
- `max_step`: maximum active-coordinate trust-region step norm.
- `prune_condition`: deterministic weak-parameter pruning target. `0.0`
  disables pruning.

`[geometry]`

- `units`: currently must be `angstrom`.
- `atoms`: ordered Cartesian atoms as `[symbol, x, y, z]`. Atom numbering is
  one-based and follows this order in constraints, substitutions and reports.
- No dummy atoms are allowed.
- No Z-matrix variables are allowed.

`[constraints]`

- `fixed_gic_patterns`: substrings matched against final working-coordinate
  labels. For `coordinate_model = "gic"` these are final GIC labels; for
  `cartesian_symmetry` they are labels such as `SC001`, `A1Cart0001` or
  `irrep=A1`. Matching parameters are reported but excluded from the fit.
- `fix_hydrogen_parameters`: optional boolean. When true, Merlino fixes a
  deterministic local coordinate frame for every H, D or T atom. The generated
  constraints use the corresponding X-H stretch, one local valence or
  linear-bend definition, and one torsional or out-of-plane orientation
  coordinate when available. This freezes hydrogen positions at the reference
  QM geometry without over-constraining heavy-atom skeletal coordinates. The
  primitive constraints are expanded over symmetry before the active GIC space
  is projected.
- `modredundant`: optional Gaussian ModRedundant freeze records. Only `F`
  records are interpreted. Supported coordinate tags are `B`, `A`, `D`, `O`
  and `L`.
- ModRedundant indexes are one-based and refer to the `[geometry].atoms` order.
  Freeze records are converted to primitive-coordinate constraints and expanded
  automatically to all symmetry-equivalent primitive coordinates before the
  active totally symmetric GIC space is projected.

`[[qm_predicates]]`

- `pattern`: substring matched against working-coordinate labels.
- `value`: target working-coordinate value in native Merlino units for GICs, or
  Angstrom Cartesian-basis amplitudes for `cartesian_symmetry`.
- `sigma`: one-sigma uncertainty in the same unit.
- `source`: free text.

`[[parameter_classes]]`

- `name`: class name reported in output tables.
- `mode`: `shared` or `fixed`.
- `patterns`: list of label substrings. With the GIC model, use explicit
  coordinate-type words such as `bond`, `angle`, `dihedral` when possible.
  With the symmetry-Cartesian model, use labels such as `SC001`, `A1Cart0001`
  or `irrep=A1`.

## 2. Gaussian Cartesian Geometry Input

Gaussian `.com/.gjf` files are accepted as geometry inputs when they contain a
Cartesian coordinate block:

```text
#p hf/sto-3g opt=modredundant

parent Cartesian geometry

0 1
C   0.000000   0.000000   0.000000
H   0.000000   0.000000   1.089000

B 1 2 F
```

Everything after the blank line following the Cartesian block is treated as
ModRedundant data. Freeze records are converted to primitive-coordinate
constraints, expanded over the detected symmetry orbit and projected onto the
active GIC space. Gaussian route sections containing `zmat` are rejected.

## 3. Isotopologue Observations

Recommended extension: `.toml`.

```toml
[[isotopologues]]
label = "parent"
substitutions = ""

[isotopologues.constants]
A_MHz = 1000.000
B_MHz = 800.000
C_MHz = 600.000

[isotopologues.vibrational_correction]
delta_A_MHz = 1.0
delta_B_MHz = 2.0
delta_C_MHz = 3.0
source = "B3LYP/cc-pVTZ"
convention = "subtract"

[isotopologues.electronic_correction]
delta_A_MHz = 0.1
delta_B_MHz = 0.2
delta_C_MHz = 0.3
source = "relativistic+BOB"
convention = "subtract"

[isotopologues.sigma_MHz]
A_MHz = 0.010
B_MHz = 0.010
C_MHz = 0.020

[[isotopologues]]
label = "13C1"
substitutions = [{ atom = 1, mass = 13 }]

[isotopologues.constants]
A_MHz = 990.000
B_MHz = 790.000
C_MHz = 590.000
```

Rules:

- `A_MHz`, `B_MHz`, `C_MHz` are experimental ground-state constants.
- `delta_*_MHz` are vibrational or electronic corrections.
- Default correction convention is `subtract`:
  `Be = B0 - Delta_vib - Delta_elec`.
- `convention = "additive"` means:
  `Be = B0 + Delta_vib + Delta_elec`.
- `substitutions` uses one-based atom indices from the parent geometry.
- `D` and `T` are accepted aliases for masses 2 and 3.
- `sigma_MHz` values are optional. If present, all three components must be
  present and positive.

CSV is supported for spreadsheet workflows. The required CSV fields are:

```text
label,A_MHz,B_MHz,C_MHz,delta_A_MHz,delta_B_MHz,delta_C_MHz,correction_source,substitutions
```

Optional CSV fields are:

```text
correction_convention,delta_elec_A_MHz,delta_elec_B_MHz,delta_elec_C_MHz,electronic_correction_source,electronic_correction_convention,sigma_A_MHz,sigma_B_MHz,sigma_C_MHz
```

## 4. Standard Outputs

Every semiexperimental run writes:

- `semiexp_report.txt`: canonical text report in `SEFIT TEXT OUTPUT v1`
  format. Required sections are `[method]`, `[constraints]`,
  `[fit_statistics]`, `[working_coordinates]`,
  `[primitive_internal_coordinates]`, `[rotational_constants]` and
  `[fit_residuals]`.
- `semiexp_report.html`: graphical self-contained report for GUI inspection.
- `semiexp_geometry.xyz`: final Cartesian equilibrium geometry.
- `semiexp_parameters.csv`: final working-coordinate parameters and propagated
  errors. These are non-redundant GICs in native GIC units for the default
  model, or Cartesian-basis amplitudes in Angstrom for the symmetry-Cartesian
  model.
- `semiexp_geometry_parameters.csv`: final bond lengths in Angstrom and angles
  or dihedrals in degrees, with propagated one-sigma errors.
- `semiexp_rotational_constants.csv`: corrected experimental constants,
  constants calculated from the fitted geometry and differences in MHz.
- `semiexp_residuals.csv`: residuals in the selected observable space.
- `semiexp_kraitchman.csv`: Kraitchman diagnostic for single substitutions.
- `semiexp_covariance.csv`, `semiexp_correlation.csv`,
  `semiexp_hessian.csv`, `semiexp_hessian_eigenvalues.csv`.
- `semiexp_diagnostics.csv`: convergence and conditioning diagnostics.
- `semiexp_influence.csv`: residual, weighted residual, chi-square contribution
  and leverage by observable.
- `semiexp_high_correlations.csv`: strongly correlated fitted-parameter pairs.
- `semiexp_manifest.json`: reproducibility manifest with input/output checksums
  and run parameters.

## 5. Command-Line Examples

Canonical job-file run:

```bash
python -m merlino semiexp \
  --job cyclopentadiene.mse.toml \
  --outdir semiexp_run
```

Interoperability run using a Gaussian Cartesian geometry file and a separate
observations file:

```bash
python -m merlino semiexp \
  --geometry cyclopentadiene_cartesian.com \
  --observations cyclopentadiene_isotopologues.toml \
  --outdir semiexp_run
```

Symmetry-Cartesian working-coordinate run:

```bash
python -m merlino semiexp \
  --job cyclopentadiene.mse.toml \
  --coordinate-model cartesian_symmetry \
  --outdir semiexp_cartesian_symmetry
```
