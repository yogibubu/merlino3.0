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

For each isotopologue the solver computes equilibrium rotational constants from
the current geometry and isotope masses. The Wilson B matrix is analytic for the
standard Merlino primitive coordinates used here (bonds, angles, linear bends,
dihedrals and out-of-plane terms). The Jacobian of rotational constants with
respect to the active non-redundant GICs is then used for the weighted
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

## Output

The output directory contains:

- `semiexp_geometry.xyz`: fitted equilibrium Cartesian geometry.
- `semiexp_parameters.csv`: final non-redundant GIC values, one-sigma errors and
  active/fixed flags.
- `semiexp_residuals.csv`: observed equilibrium constants, calculated constants
  and residuals in MHz.
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
