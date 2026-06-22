# Merlino Semiexperimental File Formats

This document defines the Merlino4 semiexperimental equilibrium-geometry file
contracts. The recommended production setup is a self-contained job file:

1. A job file, `*.mse.toml`, containing calculation keywords, Cartesian parent
   geometry, optional fixed-parameter definitions, and an inline
   `[[isotopologues]]` table with the experimental rotational data.
2. Optionally, a separate reusable isotopologue observations file, `*.toml`,
   `*.json` or `*.csv`, may be referenced from `[files].observations`.

Gaussian-style Cartesian `.com/.gjf` files remain accepted as geometry inputs
for interoperability. Legacy MSR monolithic inputs use the explicit
compatibility extensions `.msr` and `.msr.inp`; existing files named
`*_msr.inp` are accepted as a legacy alias. Generic `.inp` files are
deliberately not auto-detected as MSR files.

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

[fit]
backend = "python"
coordinate_model = "gic"
observable = "moments"
rotational_components = "auto"
max_step = 0.25
damping = 1.0e-8
step = 1.0e-4
prune_condition = 0.0
robust_loss = "none"
robust_scale = 0.0
leave_one_out = false
# optional:
# checkpoint = "semiexp_checkpoint.json"
# restart = "semiexp_checkpoint.json"

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
  "R(1,2) Frozen",
]
gic_constraints = [
  "QFIX=[GIC001+2*GIC002] Value=0.0",
  "DR(Frozen,Value=0.0)=R[1,3]-R[1,2]",
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
patterns = ["R(1,6)", "R(2,7)"]

[[parameter_classes]]
name = "XYH_angles"
mode = "fixed"
patterns = ["A("]

[[isotopologues]]
label = "parent"
[isotopologues.definition]
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
delta_A_MHz = 0.0
delta_B_MHz = 0.0
delta_C_MHz = 0.0
source = "none"
convention = "subtract"

[[isotopologues]]
label = "13C1"
[isotopologues.definition]
substitutions = [{ atom = 1, mass = 13 }]
[isotopologues.constants]
A_MHz = 990.000
B_MHz = 790.000
C_MHz = 590.000
```

### Job Tables

`[files]`

- `observations`: path to the isotopologue observations file. Relative paths
  are resolved from the job-file directory. This entry is optional when the
  job contains inline `[[isotopologues]]` tables. If both are present and no
  command-line `--observations` override is supplied, the inline tables define
  the job-local data set.

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
- `rotational_components`: `auto`, `ABC`, `AB`, `AC` or `BC`. Planar molecules
  have only two independent rotational components, so `auto` chooses one pair
  and explicit `ABC` is rejected. When `observable = "moments"`, the same
  spectroscopic labels are mapped internally to `Ia/Ib`, `Ia/Ic` or `Ib/Ic`.
  The automatic choice first minimizes the propagated instability of the
  omitted planar component and then uses Jacobian conditioning as a tie-breaker.
- `max_iter`: optional explicit iteration cap. If omitted, Merlino uses an
  automatic cap proportional to the number of effective fitted parameters.
- `step`: finite working-coordinate step used only by fallback numerical
  derivatives. Current GIC and symmetry-Cartesian SEfit paths use
  analytic Cartesian derivatives for rotational observables.
- `damping`: initial Levenberg-Marquardt damping for the trust-region solver.
- `max_step`: maximum active-coordinate trust-region step norm.
- `prune_condition`: deterministic weak-parameter pruning target. `0.0`
  disables pruning.
- `robust_loss`: optional robust IRLS loss for experimental outlier
  isotopologues. Allowed values are `none`, `huber`, `soft_l1` and `cauchy`.
  The default `none` reproduces the ordinary weighted least-squares problem.
- `robust_scale`: robust residual scale in weighted units. `0.0` selects an
  automatic median-absolute-deviation scale. Robust weights are applied only to
  experimental isotopologue rows; QM predicates keep their declared weights.
- `leave_one_out`: when `true`, Merlino performs exact leave-one-isotopologue-
  out refits after the final fit and writes `semiexp_leave_one_out.csv`.
- `checkpoint`: optional explicit checkpoint path. If omitted, runs with an
  output directory write `semiexp_checkpoint.json`.
- `restart`: optional checkpoint JSON used to restart from saved Cartesian
  coordinates. Coordinate definitions are rebuilt deterministically from the
  restart geometry.

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
- `gic_constraints`, `expression_constraints`, `fixed_expressions`,
  `gic_definitions`, `coordinate_definitions` or `definitions`: exact
  Gaussian-style expression constraints. They are appended to
  `fixed_gic_patterns` internally but are parsed as constraints, not as label
  substrings. The implemented syntax follows the Gaussian GIC keyword reference
  archived in `doc/papers/references/gaussian_gic_keywords.pdf` for the
  constraint forms needed by SEfit: `Label(Options)=Expression` and
  `Expression Options`. Accepted primitive functions are `R/B/Bond/Stretch`,
  `A/Angle/Bend`, `D/Dihedral/Torsion`, `U/out_of_plane`,
  `L/Linear/LinearBend`, `X/Y/Z`, `Cart/Cartesian` and `DotDiff`, with
  one-based atom indexes. Final GIC values can be referenced as `GIC001`,
  `GIC002`, etc. Arithmetic, parentheses, square brackets, braces, powers and
  elementary functions such as `sin`, `cos`, `arccos`, `sqrt`, `exp`, `log`,
  `min` and `max` are supported. `NAME=expression` defines a reusable
  coordinate symbol without constraining it; later expressions can reference
  that name, and cyclic definitions are rejected. `NAME(Frozen)=expression`,
  `expression Freeze` or `NAME=[expression] F` freezes the initial value of the
  expression. `Value=target` imposes an explicit Gaussian-style target in the
  SEfit constraint channel, for example `NAME(Frozen,Value=target)=expression`
  or `NAME=[expression] Value=target`.
  Targets for pure angular expressions use Gaussian default degrees when no
  unit suffix is present; this angular inference follows reusable definitions.
  `deg` and `rad` may also be given explicitly.
  Gaussian coordinate-management options such as `Inactive`, `Remove` and
  `Active` are recognized as non-constraint actions; for example
  `RPck001(Inactive,Value=...)` is not treated as a hard SEfit constraint.
  Tabulated structural constraints should always use `Value=` so that the final
  fit is not tied to the starting Cartesian geometry.
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
QFIX=[R(1,3)-R(1,2)] Value=0.0
```

Everything after the blank line following the Cartesian block is treated as
ModRedundant data. Freeze records are converted to primitive-coordinate
constraints, expanded over the detected symmetry orbit and projected onto the
active GIC space. Gaussian-style GIC expression constraints in the same block
are preserved and projected by the SEfit constraint machinery. Gaussian route
sections containing `zmat` are rejected.

## 3. Legacy MSR Compatibility Input

Standard extensions: `.msr` and `.msr.inp`. Existing files named `*_msr.inp`
are accepted as a legacy alias, but plain `.inp` is not.

This is an import format, not the recommended native Merlino format. The file
may contain either an MSR Z-matrix block followed by variable definitions, or a
direct Cartesian block with records `Atom x y z`. In both cases the geometry is
converted immediately to the canonical Merlino Cartesian representation.

For direct Cartesian MSR imports, fixed primitive parameters are supplied after
the Cartesian block with Gaussian-style ModRedundant freeze records. The block
may be introduced by `constraints`, `modredundant`, `freeze`, or `frozen`, and
may be terminated by `end`; the first isotope mass line also terminates it.
Only freeze records are interpreted:

```text
C  0.000000  0.000000  0.000000
H  0.000000  0.000000  1.089000
H  1.026719  0.000000 -0.363000

constraints
B 1 2 F
A 2 1 3 F
QFIX=[GIC001+2*GIC002] Value=0.0
end

12.000000  \parent
 1.000000
 1.000000
```

Dummy atoms (`X`, `XX`, `-1`) are allowed only in the Z-matrix branch and are
removed before GIC generation and SEfit. If a Z-matrix coordinate uses a
variable marked with `#`, for example `#RCH`, that coordinate is imported as a
frozen primitive constraint when all involved atoms are real. Coordinates
involving dummy atoms are used only to build the Cartesian geometry and are not
converted into chemical constraints.

The Z-matrix reader accepts ordinary numeric values, Fortran `D` exponents, and
variable definitions written either as `NAME = value` or `NAME value`. Atom
labels such as `C1` or `H2` are normalized to chemical symbols, while
references must always point to atoms already defined in the Z-matrix.

The isotope mass blocks define the parent and isotopologues. Labels written as
comments after a backslash, such as `\parent`, are optional; if absent, Merlino
uses deterministic labels `parent`, `iso_002`, and so on. Isotope substitutions
are inferred by comparing each mass block with the parent mass block, not from
the label text. `bexp`, `dbvib`, optional `dbelec`/`dbele`, and `weights`
sections are converted to the standard observation model. `dbvib` and
electronic corrections are interpreted with the MSR additive convention.

The same file can be used directly as a semiexperimental job:

```bash
python -m merlino semiexp --job nitrobenzene.msr.inp --outdir run_nitrobenzene
```

A complete Cartesian/ModRedundant MSR-compatible nitrobenzene example is
provided in
`doc/papers/semiexp_gic_fit/examples/nitrobenzene_cartesian_constraints.msr`.

## 4. Isotopologue Observations

Recommended location: inline in the `.mse.toml` job file. A separate `.toml`
file with the same `[[isotopologues]]` tables remains supported for reusable
data sets.

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

Inside a job file the more explicit form is preferred:

```toml
[[isotopologues]]
label = "13C1"

[isotopologues.definition]
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
- `[isotopologues.definition].substitutions` is equivalent to top-level
  `substitutions` and is preferred in self-contained job files because it
  separates the isotopologue definition from the spectroscopic data.
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
  `[constraint_diagnostics]`, `[fit_statistics]`, `[warnings]`,
  `[rank_diagnostics]`, `[working_coordinates]`,
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
- `semiexp_svd_diagnostics.csv`: singular values of the final weighted
  Jacobian and the dominant working-coordinate combinations associated with
  small singular values.
- `semiexp_uncertainty_diagnostics.csv`: covariance sensitivity to the SVD
  rank cutoff, reported as parameter sigmas for the default cutoff and several
  explicit relative cutoffs.
- `semiexp_iteration_trace.csv`: machine-readable trust-region trace for every
  accepted, rejected or topology-rejected trial step.
- `semiexp_constraints.csv`: input fixed patterns, symmetry-expanded primitive
  constraints, parameter classes and the active labels they match.
- `semiexp_warnings.csv`: non-blocking diagnostic warnings for reduced rank,
  small singular values, ill-conditioned planar component pairs, low robust
  isotopologue weights, large weighted residuals, high leverage observations,
  strongly correlated parameters, trust-region stagnation and unusually large
  propagated primitive-coordinate uncertainties.
- `semiexp_checkpoint.json`: restartable solver state with Cartesian geometry,
  damping/trust radius, labels, active labels and robust weights.
- `semiexp_leave_one_out.csv`: exact leave-one-isotopologue-out refits, written
  only when `leave_one_out = true`.
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
