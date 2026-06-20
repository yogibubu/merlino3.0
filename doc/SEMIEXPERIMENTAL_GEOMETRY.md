# Semiexperimental Equilibrium Geometries

This is the Merlino4 standard solver for semiexperimental equilibrium
geometries. It deliberately avoids hand-built structural parameterizations:
the input geometry is Cartesian, the optimized
parameters are non-redundant Merlino GICs, and the final result is a Cartesian
equilibrium structure with propagated errors for the fitted internal
parameters.

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
- Analytic Wilson B matrix for standard internal primitives.
- Weighted Levenberg-Marquardt least squares with adaptive damping.
- Direct propagation of experimental uncertainties to GIC parameters.
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
- `--max-step 0.25` limits the norm of active-GIC steps and prevents aggressive
  updates from leaving the chemically valid topology basin.
- `--prune-condition 0` leaves the full totally symmetric GIC subspace active.
  Passing a positive value enables deterministic removal of weak A1 parameters
  until the initial weighted Jacobian condition is below the requested target.
  Removed coordinates are reported as `auto_pruned_weak`; this is a diagnostic
  fallback, not the recommended default for MSR-grade refinements.
- `--damping 1e-8` is only the initial Levenberg-Marquardt damping. It is
  decreased after accepted steps and increased after rejected steps.
- `--max-iter` defaults to an automatic cap proportional to the number of
  effective optimized parameters: `max(8, 2*N)`. Passing a positive value keeps
  an explicit user cap.

Use `--observable rotational_constants` only when the scientific comparison
must be made directly in MHz.

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
`B/A/D/O/L ... F` are converted into fixed GIC patterns before the fit.

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
- moment or rotational-constant target;
- rotational-constant component policy;
- fixed GIC patterns;
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
9. Compute the Jacobian of observables with respect to active GICs.
10. Solve weighted LM normal equations with adaptive damping and step limiting.
11. Back-transform GIC steps to Cartesian displacements using the analytic B
   matrix and reject steps that do not improve the weighted objective.
12. Recompute covariance, correlation, Hessian eigenvalues and diagnostics at
   the final geometry.

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
least-squares parameters through the same GIC-to-Cartesian differential used by
the optimizer:

```text
J_top = B_top pinv(B_GIC) T_active
Cov(topological) = J_top Cov(active) J_top^T
```

`B_top` is the Wilson B matrix of the reported bond, angle and dihedral
coordinates, `B_GIC` is the B matrix of the final non-redundant GIC set and
`T_active` is the active parameter/class transform. Bond sigmas are reported in
Angstrom; angle and dihedral sigmas are reported in degrees.

## Output

The output directory contains:

- `semiexp_geometry.xyz`: fitted equilibrium Cartesian geometry.
- `semiexp_report.txt`: plain-text human-readable report with diagnostics,
  corrected/calculated rotational constants, final bond lengths/angles/dihedrals
  with propagated errors, GIC parameters and residuals.
- `semiexp_report.html`: self-contained run report with diagnostics, parameter
  classes, fitted GICs, final Cartesian bond lengths/angles/dihedrals with
  propagated errors, rotational-constant comparison, residuals and Kraitchman
  comparison.
- `semiexp_tables.tex`: paper-ready LaTeX tabular fragments for parameters,
  rotational constants, residuals and Kraitchman comparison.
- `semiexp_parameters.csv`: final non-redundant GIC values, one-sigma errors and
  active/fixed flags.
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

The parameter values use native Merlino GIC units: stretches in Angstrom and
angular coordinates in radians.

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
