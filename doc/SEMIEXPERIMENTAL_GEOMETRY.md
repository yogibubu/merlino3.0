# Semiexperimental Equilibrium Geometries

Merlino4 fits semiexperimental equilibrium geometries from a parent Cartesian
structure and isotopologue rotational constants. The solver is independent from
Gaussian: Gaussian can be used upstream to provide vibrational corrections, but
the fit consumes only Merlino data files.

## Input

Run from the CLI with:

```bash
python -m merlino semiexp \
  --xyz parent_initial.xyz \
  --observations isotopologues.csv \
  --outdir semiexp_run \
  --observable moments \
  --max-step 0.25
```

The XYZ file contains the starting parent geometry in Angstrom.

The observation CSV columns are:

```text
label,A_MHz,B_MHz,C_MHz,delta_A_MHz,delta_B_MHz,delta_C_MHz,correction_source,substitutions
```

`A_MHz`, `B_MHz` and `C_MHz` are experimental ground-state constants `B0`.
`delta_*_MHz` are vibrational corrections in the same convention used by
Merlino:

```text
Be = B0 - delta
```

`substitutions` is a semicolon-separated list of one-based atom substitutions,
for example `2:13;5:18`. Empty substitutions mean the parent isotopologue.

Optional columns `sigma_A_MHz`, `sigma_B_MHz` and `sigma_C_MHz` provide
experimental uncertainties. When present, Merlino uses inverse-variance weights
`1/sigma^2` and propagates these uncertainties to the fitted GIC parameters.

## Fit Model

Merlino generates primitive internal coordinates from the starting Cartesian
geometry, builds the same non-redundant GIC transform used by the GF workflow,
and optimizes active GIC values by least squares.

For each isotopologue the solver computes either principal moments of inertia
or rotational constants from the current geometry and isotope masses. The
default is `--observable moments`, because moments are linear in mass geometry
and avoid the reciprocal amplification present in rotational constants. Use
`--observable rotational_constants` only when that is the intended experimental
fit target.

For planar molecules and `--observable rotational_constants
--rotational-components auto`, Merlino evaluates the initial Jacobian for
`AB`, `AC` and `BC` and selects the pair with best rank and smallest condition
number. Non-planar molecules use `ABC` by default.

The Wilson B matrix is analytic for the standard Merlino primitive coordinates
used here (bonds, angles, linear bends, dihedrals and out-of-plane terms). The
Jacobian with respect to active non-redundant GICs is used for the weighted
least-squares normal equations and for error propagation.

Parameters can be frozen with:

```bash
python -m merlino semiexp ... --fixed "GIC001,angle"
```

Each token is matched as a case-insensitive substring of the generated GIC
labels. Fixed parameters are reported but excluded from the least-squares
normal equations.

The optimizer uses a Levenberg-Marquardt style weighted least-squares step with
adaptive damping. Steps that do not improve the weighted objective are rejected,
the damping is increased, and the next iteration retries a more conservative
normal equation. `--max-step` limits the active-GIC step norm and is useful when
the starting geometry is only approximate.

QM-estimated parameters can be added as weighted predicates:

```bash
python -m merlino semiexp ... --qm-predicate "GIC001:1.234:0.010:qm"
```

The format is `label_pattern:value:sigma[:source]`. The label pattern is matched
against generated GIC labels and contributes a pseudo-observation with weight
`1/sigma^2`.

## Output

The output directory contains:

- `semiexp_geometry.xyz`: fitted equilibrium Cartesian geometry.
- `semiexp_parameters.csv`: final non-redundant GIC values, one-sigma errors and
  active/fixed flags.
- `semiexp_residuals.csv`: observed, calculated and residual values for the
  selected observable. Units are MHz for rotational constants, amu Angstrom^2
  for moments of inertia and native GIC units for QM predicates.
- `semiexp_covariance.csv`: propagated covariance matrix for active parameters.
- `semiexp_correlation.csv`: correlation matrix for active parameters.
- `semiexp_hessian.csv`: Gauss-Newton least-squares Hessian.
- `semiexp_hessian_eigenvalues.csv`: eigenvalues used to classify the fitted
  stationary point as `minimum`, `flat_or_rank_deficient` or
  `transition_state_or_saddle`.
- `semiexp_diagnostics.csv`: convergence reason, objective, weighted RMS,
  reduced chi square, Jacobian rank, condition number and accepted/rejected
  steps.
- `semiexp_manifest.json`: reproducibility manifest with checksums.

The parameter values use the native Merlino GIC units: stretches in Angstrom and
angular coordinates in radians.
