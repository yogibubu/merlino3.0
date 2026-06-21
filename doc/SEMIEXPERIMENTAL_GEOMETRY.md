# Semiexperimental Equilibrium Geometries

This is the Merlino4 standard solver for semiexperimental equilibrium
geometries. It deliberately avoids hand-built structural parameterizations:
the input geometry is Cartesian, the optimized
parameters are non-redundant Merlino GICs by default, and the final result is
a Cartesian equilibrium structure with propagated errors for the fitted
internal parameters. The alternative working-coordinate model uses
Hessian-free symmetry-adapted Cartesian displacements; the final structural
report is still given as primitive bond lengths, angles and dihedrals with
propagated errors.

## Why This Is The Standard Solver

Classical semiexperimental geometry programs usually require a hand-built
internal-coordinate template. That is fragile because the result depends on
coordinate ordering, dummy atoms, manually chosen dependent coordinates and
molecule-specific parameter choices. The Merlino solver instead uses:

- Cartesian parent geometry as the only structural input.
- Automatic topology and primitive GIC generation.
- Reduction of the generated GICs to a non-redundant active set.
- Automatic point-group identification and symmetry adaptation of the
  non-redundant GICs within homogeneous coordinate families.
- Analytic Wilson B matrix for standard internal primitives and analytic
  Cartesian derivatives of principal moments/rotational constants.
- Weighted trust-region Levenberg-Marquardt least squares with coordinate-block
  scaling and predicted/actual reduction control.
- Direct propagation of experimental uncertainties to GIC parameters.
- Optional totally symmetric symmetry-Cartesian working coordinates requiring
  neither a Hessian nor a Wilson B-matrix inversion.
- Optional QM predicates as weighted priors, not hard constraints.
- Manifested outputs and diagnostics suitable for regression checks.

This makes the workflow more general for rings, fused systems, bridge atoms,
planar molecules and cases where a conventional hand-built coordinate template
would be ambiguous or ill-conditioned.

## Recommended Defaults

The recommended CLI is:

```bash
python -m merlino semiexp \
  --job cyclopentadiene.mse.toml \
  --outdir semiexp_run
```

For reduced-dimensionality fits without deuterium substitutions, the local
geometry of every H, D or T atom can be fixed in one operation:

```bash
python -m merlino semiexp \
  --job norcamphor.mse.toml \
  --fix-hydrogens \
  --outdir semiexp_norcamphor
```

The equivalent job-file setting is `fix_hydrogen_parameters = true` under
`[constraints]`.

The equivalent interoperability mode, useful when the parent geometry is
already available as a Gaussian Cartesian input, is:

```bash
python -m merlino semiexp \
  --geometry parent_initial.com \
  --observations isotopologues.toml \
  --outdir semiexp_run \
  --observable moments \
  --rotational-components auto \
  --max-step 0.25
```

These defaults are intentional:

- `--observable moments` is the standard target. Principal moments of inertia
  are more stable than rotational constants because rotational constants are
  reciprocal in the moments and amplify small errors when moments are small.
- `--rotational-components auto` only matters if `--observable
  rotational_constants` is selected. For planar molecules it selects the
  best-conditioned pair among `AB`, `AC` and `BC`; non-planar molecules use
  `ABC`.
- `--max-step 0.25` limits the norm of active-coordinate steps and prevents
  aggressive updates from leaving the chemically valid topology basin.
- `--prune-condition 0` leaves the full totally symmetric GIC subspace active.
  Passing a positive value enables deterministic removal of weak A1 parameters
  until the initial weighted Jacobian condition is below the requested target.
  Removed coordinates are reported as `auto_pruned_weak`; this is a diagnostic
  fallback, not the recommended default for MSR-grade refinements.
- `--damping 1e-8` is only the initial Levenberg-Marquardt damping. It is
  updated from the ratio between predicted and actual objective reduction, not
  from a blind accept/reject multiplier.
- `--max-iter` defaults to an automatic cap proportional to the number of
  effective optimized parameters: `max(8, 2*N)`. Passing a positive value keeps
  an explicit user cap.

Use `--observable rotational_constants` only when the scientific comparison
must be made directly in MHz.

The default coordinate model is:

```bash
python -m merlino semiexp --job parent.mse.toml --coordinate-model gic
```

The coordinate model is a scientific choice, not only an implementation detail.
For experimental structure refinement the ideal working coordinates are
non-redundant, symmetry-adapted, compatible with constraints and parameter
classes, and suitable for rigorous covariance propagation to ordinary
structural parameters. Rotational spectroscopy is the current Merlino working
example, but the same requirements apply to least-squares refinement against
any experimental structural observable.

The symmetry-Cartesian model is the B-free alternative:

```bash
python -m merlino semiexp \
  --job parent.mse.toml \
  --coordinate-model cartesian_symmetry \
  --outdir semiexp_cartesian_symmetry
```

This model uses only the parent Cartesian geometry. Merlino projects out
translations and rotations, applies point-group symmetry projectors to the
remaining Cartesian displacement space, and optimizes only the totally
symmetric directions. It does not require a Hessian, a force-field calculation
or a Wilson B-matrix pseudoinverse. Primitive constraints are still enforced
through analytic primitive derivatives, and the final errors are propagated to
ordinary internal coordinates from the fitted Cartesian covariance.

## Input

The canonical Merlino input is a two-file setup:

- `*.mse.toml`: job file with keywords, Cartesian parent geometry and optional
  fixed-parameter definitions;
- `*.toml`, `*.json` or `*.csv`: isotopologue observations with rotational
  constants, vibrational corrections, electronic corrections and optional
  uncertainties.

The complete file-format contract is documented in
`doc/SEMIEXPERIMENTAL_FILE_FORMATS.md`.

For interoperability, the parent geometry can still be provided as standard XYZ
or as a Gaussian `.com` / `.gjf` input with Cartesian coordinates. Gaussian
Z-matrices are intentionally not accepted in this workflow. If the Gaussian
input contains a ModRedundant section, freeze constraints written as
`B/A/D/O/L ... F` are converted into primitive-coordinate constraints. Each
fixed primitive is expanded automatically to all symmetry-equivalent primitives
before the constraint projector is built in the active GIC space.
For limited-isotopologue data sets, `[constraints] fix_hydrogen_parameters =
true` builds a deterministic local coordinate frame for every H, D or T atom
and fixes those primitive constraints in the same projected-constraint
machinery. This keeps hydrogen positions at their QM reference values without
using the full over-complete set of hydrogen-containing primitives.

Minimal accepted Gaussian-style geometry input:

```text
#p hf/sto-3g opt=modredundant

parent Cartesian geometry

0 1
C   0.000000   0.000000   0.000000
H   0.000000   0.000000   1.089000

B 1 2 F
```

Everything after the blank line following the Cartesian block is interpreted as
Gaussian ModRedundant data. Fixed parameters must be expressed there with the
standard freeze action `F`; Z-matrix variables and MSR-style `R0001/A0001`
cards are not part of the Merlino4 SE input.

Isotopologue observations are provided separately as TOML, JSON or CSV. TOML is
the recommended human-edited format because it keeps constants, isotope
substitutions and corrections grouped by isotopologue.

Recommended TOML:

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
label = "D2"
substitutions = [{ atom = 2, mass = 2 }]

[isotopologues.constants]
A_MHz = 900.000
B_MHz = 700.000
C_MHz = 500.000
```

The same schema can be written as JSON for GUI/script generation:

```json
{
  "isotopologues": [
    {
      "label": "13C1",
      "substitutions": {"1": 13},
      "constants": {"A_MHz": 990.0, "B_MHz": 790.0, "C_MHz": 590.0},
      "vibrational_correction": {
        "delta_A_MHz": 0.5,
        "delta_B_MHz": 1.5,
        "delta_C_MHz": 2.5,
        "source": "QM"
      },
      "electronic_correction": {
        "delta_A_MHz": 0.0,
        "delta_B_MHz": 0.0,
        "delta_C_MHz": 0.0,
        "source": "none"
      }
    }
  ]
}
```

CSV remains supported for spreadsheets and backward compatibility.

The observation CSV columns are:

```text
label,A_MHz,B_MHz,C_MHz,delta_A_MHz,delta_B_MHz,delta_C_MHz,correction_source,correction_convention,substitutions
```

`A_MHz`, `B_MHz` and `C_MHz` are experimental ground-state constants `B0`.
`delta_*_MHz` are vibrational corrections. Optional electronic corrections are
available in TOML/JSON as `electronic_correction` and in CSV as
`delta_elec_A_MHz`, `delta_elec_B_MHz`, `delta_elec_C_MHz`,
`electronic_correction_source` and `electronic_correction_convention`.

Merlino defaults to the subtractive convention:

```text
Be = B0 - delta_vib - delta_elec
```

Set `convention = "additive"` for MSR-style corrections already defined as
terms to add to the ground-state constants:

```text
Be = B0 + delta_vib + delta_elec
```

`substitutions` uses one-based atom indices. In TOML/JSON it can be a mapping,
a list of `{atom, mass}` records, or the compact string used by CSV. `D` and
`T` are accepted aliases for masses 2 and 3. Empty substitutions mean the
parent isotopologue.

Optional columns `sigma_A_MHz`, `sigma_B_MHz` and `sigma_C_MHz` provide
experimental uncertainties. When present, Merlino uses inverse-variance weights
`1/sigma^2`. If the fit target is moments, these uncertainties are propagated
through `I = K/B`.

## GUI Workflow

The Merlino4 dashboard exposes the semiexperimental solver from the
`Semiexperimental Geometry` workflow. The panel lets the user select:

- parent Cartesian geometry file (`.xyz`, `.com` or `.gjf`);
- TOML/JSON/CSV isotopologue observations;
- or a complete Merlino job file (`.mse.toml`);
- output directory;
- Python or Fortran77 backend request;
- GIC or symmetry-Cartesian coordinate model;
- moment or rotational-constant target;
- rotational-constant component policy;
- fixed GIC patterns;
- automatic fixing of local hydrogen/deuterium/tritium geometry constraints;
- QM predicate observations;
- shared or fixed parameter classes.

The GUI builds the same command used by the CLI and runs it asynchronously, so
the interface remains responsive during the least-squares fit.

The same panel also provides operational helpers:

- an isotopologue table editor that writes the recommended TOML input;
- a structured GIC preview showing labels, coordinate type, atom indices,
  automatic class suggestions and active/fixed state before the fit;
- input validation for duplicate labels, impossible substitutions, suspicious
  corrections and class definitions that match no GIC or mix coordinate types;
- conditioning preview for the selected isotopologues, fixed parameters and
  classes before launching the least-squares fit;
- automatic suggestions for shared/fixed parameter classes;
- save/load of GUI presets for reproducible setup of repeated fits;
- direct opening of the generated `semiexp_report.html`.
- guided workflow state from manifests and expected output files.

## Fit Model

The default GIC fit model is:

1. Read the parent Cartesian geometry (`.xyz`, `.com` or `.gjf`).
2. Build topology and primitive internal coordinates.
3. Build the primitive GIC set automatically, as in other black-box internal
   coordinate generators.
4. Keep every stretching coordinate as its primitive bond `R(i,j)`. Reduce
   bending, linear-bending, torsional, out-of-plane and ring candidate sets to
   a non-redundant transform.
5. Symmetrize and label the resulting coordinates after automatic point-group
   identification.
6. For rings, use endocyclic dihedral combinations and analogous endocyclic
   valence-angle combinations; Cremer-Pople variables are not used as GIC fit
   coordinates.
7. Convert the selected observations to the fit target:
   moments of inertia by default, rotational constants on request.
8. Add optional QM predicates as weighted pseudo-observations.
9. Compute the Jacobian of observables with respect to active working
   coordinates from analytic Cartesian derivatives of principal moments or
   rotational constants. The finite-difference path is retained only as a
   parallel fallback for future non-analytic observables.
10. Solve weighted trust-region LM equations with Cauchy fallback, predicted
    objective reduction, homogeneous coordinate-block scaling and adaptive
    trust-radius/damping control.
11. For the GIC model, back-transform GIC steps to Cartesian displacements
    using the analytic B matrix. Line-search trials reuse the current GICForge
    coordinate model. Every accepted GIC step is validated by rerunning
    GICForge and comparing the point-group, irrep and coordinate-family
    signature with the reference model; topology-changing steps are rejected
    and the trust radius is reduced. The Cartesian-GIC projector is refreshed
    analytically only when needed and
    otherwise updated by a secant correction.
12. Recompute covariance, correlation, Hessian eigenvalues and diagnostics at
    the final geometry.

The symmetry-Cartesian fit model shares the same observation model,
trust-region least-squares solver, parameter classes, primitive constraints,
covariance analysis and final reporting. Its coordinate stage is:

1. Orient the parent Cartesian geometry and detect the point group.
2. Build the Cartesian vibrational projector by removing translations and
   rotations.
3. Project the remaining displacement space by point-group irreps.
4. Keep only totally symmetric Cartesian displacements as active fit
   coordinates.
5. Update Cartesian geometry directly as linear combinations of the selected
   symmetry-Cartesian vectors; no GIC B projector is required.
6. Evaluate final topological bond lengths, angles and dihedrals from the
   optimized Cartesian structure and propagate their errors through the
   Cartesian-basis differential.

The run directory includes `semiexp_influence.csv` for residual, weighted
residual, chi-square contribution and leverage of every observable, and
`semiexp_high_correlations.csv` for strongly correlated fitted parameters. The
GUI prints these paths and the solver counters in its expert diagnostics block
after a run.

`scripts/run_semiexp_benchmarks.py` runs the standard regression benchmark set:
water, cyclopentadiene, azulene, uracil, pyridine and guanine for GIC/symmetry
checks, plus the semiexperimental fit cases that have observations in the repo.
The azulene fit benchmark reads the MSR-like `.mse.toml` job so that the same
primitive ModRedundant constraints and numerical settings used in production are
exercised.

The Wilson B matrix is analytic for Merlino's standard primitives: bonds,
angles, linear bends, dihedrals and out-of-plane terms. Fragment coordinates
retain their existing finite-difference fallback, but semiexperimental
molecular GIC fits use connected molecular coordinates.

The non-redundant and symmetrized GICs are built block-wise. Stretches are
direct primitive bond coordinates. Valence angles, linear bends, torsions,
out-of-plane terms and ring-specific combinations are pruned and symmetrized
within their own coordinate families. This prevents the least-squares variables
from mixing physically different coordinate types.

## QM Predicates

QM-estimated parameters can be included as weighted priors:

```bash
python -m merlino semiexp ... --qm-predicate "GIC001:1.234:0.010:qm"
```

The format is:

```text
label_pattern:value:sigma[:source]
```

The label pattern is matched against generated GIC labels. Each match adds a
pseudo-observation with weight `1/sigma^2`. Predicates are soft constraints:
they stabilize underdetermined fits without hiding disagreement between
experiment and the QM estimate.

## Fixed Parameters And Parameter Classes

Parameters can be frozen with:

```bash
python -m merlino semiexp ... --fixed "GIC001,angle"
```

Each token is matched as a case-insensitive substring of the generated GIC
labels. Fixed parameters are reported but excluded from the least-squares
normal equations.

Classes of geometric parameters can also be corrected together or blocked
together:

```bash
python -m merlino semiexp ... \
  --parameter-class "CH:shared:bond(1,2)|bond(1,3)" \
  --parameter-class "XYH:fixed:angle"
```

The format is:

```text
name:shared|fixed:pattern[|pattern...]
```

`shared` compresses all matched active GICs to one least-squares variable and
applies the same correction to the whole class. This is useful for CH stretches
or chemically equivalent XH distances when the corresponding deuterated
isotopologues are missing. `fixed` blocks the complete class, which is useful
for XYH or XH2 angle families that should be retained from the starting/QM
model. Class matching is type-safe in practice because patterns are matched on
explicit generated GIC labels; users should include enough label context
(`bond`, `angle`, atom indices or family-specific substrings) to avoid mixing
different coordinate types.

The covariance and Hessian are computed for the effective least-squares
variables. The reported parameter table expands the resulting one-sigma error
back to each member of a shared class and records the class name.

## Kraitchman Comparison

For isotopologues with exactly one substitution, Merlino writes
`semiexp_kraitchman.csv`. The file compares absolute substitution coordinates
from the Kraitchman equations with the absolute coordinates of the fitted
structure in the parent principal-axis frame.

The implementation is in `merlino_semiexp.kraitchman`. It uses corrected
equilibrium rotational constants, converts them to moments of inertia, and
uses the standard reduced substitution mass
`Delta m * M_parent / (M_parent + Delta m)`. The sign of each coordinate is
assigned from the current fitted geometry in the parent principal-axis frame.
When single-substitution data are available, Merlino also writes
`semiexp_kraitchman_geometry.xyz`, a Kraitchman-seeded geometry in the same
principal-axis frame. With at least three well-conditioned substituted atoms,
the seed applies a rigid Kabsch update before overwriting the substituted atom
positions; otherwise it directly updates only the substituted atoms.

This comparison is diagnostic, not a replacement for the semiexperimental fit:
Kraitchman coordinates lose signs, are most informative for single substitutions
and do not exploit the full correlated least-squares model. They are useful for
spotting inconsistent assignments, problematic vibrational corrections or
outlier isotopologues before accepting the final fit.

## Planar Molecules

For planar molecules, fitting all three rotational constants can be less stable
than fitting a well-conditioned pair. With:

```bash
--observable rotational_constants --rotational-components auto
```

Merlino evaluates the initial Jacobian for `AB`, `AC` and `BC`, then chooses
the pair with highest rank and lowest condition number. The selected components
are written to `semiexp_diagnostics.csv`.

The default `--observable moments` remains preferred for planar and non-planar
molecules unless direct MHz residuals are specifically required.

## Weak-Parameter Pruning

GICForge generates a complete totally symmetric non-redundant coordinate set,
but a given isotopologue set does not necessarily observe every A1 direction
with comparable accuracy. Leaving a weak direction active can inflate the
standard deviations of the fitted coordinates without improving the residual.

Before the nonlinear fit starts, Merlino analyzes the weighted Jacobian of the
selected observables with respect to the active A1 coordinates. If the condition
number is above `--prune-condition`, it removes the coordinate column whose
removal most improves the condition number, repeats deterministically, and then
uses the pruned active set for the fit. The default target is 200. Pruned
coordinates remain in `semiexp_parameters.csv` with `active=0` and
`parameter_class=auto_pruned_weak`; diagnostics and the manifest record the
exact patterns removed.

## Topological Parameter Uncertainties

The final human-readable structural table is not a second fit. Merlino evaluates
ordinary topological coordinates on the optimized Cartesian geometry: bonded
distances, valence angles and proper dihedrals from the final connectivity.
Their one-sigma errors are propagated from the covariance matrix of the active
least-squares parameters through the same working-coordinate-to-Cartesian
differential used by the optimizer. For the GIC model this is:

```text
J_top = B_top pinv(B_GIC) T_active
Cov(topological) = J_top Cov(active) J_top^T
```

`B_top` is the Wilson B matrix of the reported bond, angle and dihedral
coordinates, `B_GIC` is the B matrix of the final non-redundant GIC set and
`T_active` is the active parameter/class transform. Bond sigmas are reported in
Angstrom; angle and dihedral sigmas are reported in degrees.

For the symmetry-Cartesian model, `pinv(B_GIC) T_active` is replaced by the
symmetry-Cartesian basis times the active parameter/class transform. The
reported primitive internal coordinates and their errors therefore have the
same interpretation in both coordinate models.

## Output

The output directory contains:

- `semiexp_geometry.xyz`: fitted equilibrium Cartesian geometry.
- `semiexp_report.txt`: canonical plain-text `SEFIT TEXT OUTPUT v1` report. It
  records method, solver, coordinate model and coordinate basis used in the fit,
  input constraints, QM predicates and parameter classes, fit statistics,
  working-coordinate values/errors, primitive internal bond lengths, angles and
  dihedrals with propagated errors, rotational-constant comparison and fit
  residuals.
- `semiexp_report.html`: self-contained run report with diagnostics, parameter
  classes, fitted GICs, final Cartesian bond lengths/angles/dihedrals with
  propagated errors, rotational-constant comparison, residuals and Kraitchman
  comparison.
- `semiexp_tables.tex`: paper-ready LaTeX tabular fragments for parameters,
  rotational constants, residuals and Kraitchman comparison.
- `semiexp_parameters.csv`: final working-coordinate values, one-sigma errors
  and active/fixed flags. These are non-redundant GIC values for
  `coordinate_model = "gic"` and Cartesian-basis amplitudes for
  `coordinate_model = "cartesian_symmetry"`.
- `semiexp_geometry_parameters.csv`: final Cartesian geometry interpreted as
  ordinary structural parameters. Bond lengths are reported in Angstrom;
  valence angles and proper dihedrals are reported in degrees. The table also
  contains propagated one-sigma errors from the semiexperimental covariance
  matrix, atom indices and element labels.
- `semiexp_residuals.csv`: observed, calculated and residual values for the
  selected observable. Units are MHz for rotational constants, amu Angstrom^2
  for moments of inertia and native GIC units for QM predicates.
- `semiexp_rotational_constants.csv`: for every isotopologue and for A, B and C,
  the corrected experimental equilibrium rotational constants, the constants
  calculated from the fitted geometry and the difference
  `corrected experimental - calculated`, all in MHz. This file is always
  written, even when the fit target is moments of inertia.
- `semiexp_kraitchman.csv`: diagnostic comparison with Kraitchman substitution
  coordinates for single-substitution isotopologues.
- `semiexp_kraitchman_geometry.xyz`: diagnostic Kraitchman-seeded geometry in
  the parent principal-axis frame, written when single-substitution
  isotopologues are available.
- `semiexp_covariance.csv`: propagated covariance matrix for active parameters.
- `semiexp_correlation.csv`: correlation matrix for active parameters.
- `semiexp_hessian.csv`: Gauss-Newton least-squares Hessian.
- `semiexp_hessian_eigenvalues.csv`: eigenvalues used to classify the fitted
  stationary point as `minimum`, `flat_or_rank_deficient` or
  `transition_state_or_saddle`.
- `semiexp_diagnostics.csv`: convergence reason, objective, weighted RMS,
  reduced chi square, Jacobian rank, condition number, accepted/rejected steps,
  automatic/explicit iteration cap, selected observable and selected
  components, plus any `auto_pruned_weak` parameters removed from the active
  fit.
- `semiexp_manifest.json`: reproducibility manifest with checksums.

The GIC parameter values use native Merlino units: stretches in Angstrom and
angular coordinates in radians. Normal-mode working parameters are Cartesian
displacement amplitudes in Angstrom along normalized mode vectors.

The manifest records more than file paths: backend role, Fortran77 kernel
source, GIC generation policy, isotopologue count, predicate count, active and
effective parameter counts, iteration cap, rank, condition number, weighted
RMS, reduced chi-square, parameter classes and ring-coordinate convention.

## Benchmarks

Benchmark runs are represented by `SemiexperimentalBenchmarkCase` objects and
summarized with `run_semiexperimental_benchmark`. The resulting CSV table
records RMS, iteration count, rank, condition number, stationary-point
classification, number of fitted parameters and available Kraitchman rows. This
is the recommended format for MSR-style validation sets and ring/fused-ring
stress tests.

The repository also contains `benchmarks/semiexp_msr/manifest.toml`, which is
the schema for curated MSR-style validation cases. Each case should record the
parent Cartesian geometry, the isotopologue observation table, the reference source and the
expected numerical diagnostics.

Small executable examples are under `examples/semiexp/`:

```bash
python -m merlino semiexp \
  --xyz examples/semiexp/water/parent.xyz \
  --observations examples/semiexp/water/isotopologues.toml \
  --outdir working/examples/water_semiexp
```

The same files can be selected from the GUI semiexperimental workflow.

## Fortran77 Role And Merlino3 Regression

The semiexperimental production workflow is Python-orchestrated. The Fortran77
semiexp source is a validated numerical-kernel layer for analytic B rows,
rotational constants and least-squares normal equations, including compressed
parameter classes. It is deliberately not a second metadata/input/reporting
implementation.

GIC regression against the frozen Merlino3 baseline is handled by:

```bash
python scripts/compare_gic_merlino3.py --fixture path/to/gic_fixture --out gic_regression.json
```

The fixture must contain `provin`. The script runs Merlino3 and Merlino4
GICForge executables when available, normalizes selected text outputs and
reports exact matches/mismatches. Without `--strict`, missing local baselines
are reported as skipped rather than failing routine CI.

Granular freeze targets are available for local checks:

```bash
scripts/freeze_gui.sh
scripts/freeze_fortran.sh
scripts/freeze_semiexp.sh
scripts/freeze_regression.sh
```

## Quality Checks

A production fit should be accepted only after checking:

- `stationary_point` is `minimum`.
- Jacobian rank is sufficient for the number of active parameters.
- Condition number is not pathologically large.
- Correlations do not show near-linear parameter dependence unless chemically
  expected.
- Residuals are compatible with experimental uncertainties and QM vibrational
  correction quality.
- Any QM predicates have residuals consistent with their assigned sigma.

If the fit is rank-deficient, prefer adding more isotopologues, reducing active
parameters, or adding physically justified QM predicates. Do not use hard
constraints just to hide an underdetermined model.
