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
- Automatic topology and non-redundant GIC generation.
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
  --xyz parent_initial.xyz \
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
- `--damping 1e-8` is only the initial Levenberg-Marquardt damping. It is
  decreased after accepted steps and increased after rejected steps.

Use `--observable rotational_constants` only when the scientific comparison
must be made directly in MHz.

## Input

The XYZ file contains the starting parent geometry in Angstrom. Isotopologue
observations can be provided as TOML, JSON or CSV. TOML is the recommended
human-edited format because it keeps constants, isotope substitutions and
corrections grouped by isotopologue.

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

[isotopologues.electronic_correction]
delta_A_MHz = 0.1
delta_B_MHz = 0.2
delta_C_MHz = 0.3
source = "relativistic+BOB"

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
label,A_MHz,B_MHz,C_MHz,delta_A_MHz,delta_B_MHz,delta_C_MHz,correction_source,substitutions
```

`A_MHz`, `B_MHz` and `C_MHz` are experimental ground-state constants `B0`.
`delta_*_MHz` are vibrational corrections. Optional electronic corrections are
available in TOML/JSON as `electronic_correction` and in CSV as
`delta_elec_A_MHz`, `delta_elec_B_MHz`, `delta_elec_C_MHz`,
`electronic_correction_source`.

Merlino uses the convention:

```text
Be = B0 - delta_vib - delta_elec
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

- parent Cartesian XYZ file;
- TOML/JSON/CSV isotopologue observations;
- output directory;
- Python or Fortran77 backend request;
- moment or rotational-constant target;
- rotational-constant component policy;
- fixed GIC patterns;
- QM predicate observations;
- shared or fixed parameter classes.

The GUI builds the same command used by the CLI and runs it asynchronously, so
the interface remains responsive during the least-squares fit.

## Fit Model

1. Read the parent XYZ geometry.
2. Build topology and primitive internal coordinates.
3. Build the non-redundant GIC transform used by Merlino GF workflows.
4. Convert the selected observations to the fit target:
   moments of inertia by default, rotational constants on request.
5. Add optional QM predicates as weighted pseudo-observations.
6. Compute the Jacobian of observables with respect to active GICs.
7. Solve weighted LM normal equations with adaptive damping and step limiting.
8. Back-transform GIC steps to Cartesian displacements using the analytic B
   matrix and reject steps that do not improve the weighted objective.
9. Recompute covariance, correlation, Hessian eigenvalues and diagnostics at
   the final geometry.

The Wilson B matrix is analytic for Merlino's standard primitives: bonds,
angles, linear bends, dihedrals and out-of-plane terms. Fragment coordinates
retain their existing finite-difference fallback, but semiexperimental
molecular GIC fits use connected molecular coordinates.

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

## Output

The output directory contains:

- `semiexp_geometry.xyz`: fitted equilibrium Cartesian geometry.
- `semiexp_parameters.csv`: final non-redundant GIC values, one-sigma errors and
  active/fixed flags.
- `semiexp_residuals.csv`: observed, calculated and residual values for the
  selected observable. Units are MHz for rotational constants, amu Angstrom^2
  for moments of inertia and native GIC units for QM predicates.
- `semiexp_kraitchman.csv`: diagnostic comparison with Kraitchman substitution
  coordinates for single-substitution isotopologues.
- `semiexp_covariance.csv`: propagated covariance matrix for active parameters.
- `semiexp_correlation.csv`: correlation matrix for active parameters.
- `semiexp_hessian.csv`: Gauss-Newton least-squares Hessian.
- `semiexp_hessian_eigenvalues.csv`: eigenvalues used to classify the fitted
  stationary point as `minimum`, `flat_or_rank_deficient` or
  `transition_state_or_saddle`.
- `semiexp_diagnostics.csv`: convergence reason, objective, weighted RMS,
  reduced chi square, Jacobian rank, condition number, accepted/rejected steps,
  selected observable and selected components.
- `semiexp_manifest.json`: reproducibility manifest with checksums.

The parameter values use native Merlino GIC units: stretches in Angstrom and
angular coordinates in radians.

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
