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
  --outdir semiexp_run
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

## Fit Model

Merlino generates primitive internal coordinates from the starting Cartesian
geometry, builds the same non-redundant GIC transform used by the GF workflow,
and optimizes active GIC values by least squares.

For each isotopologue the solver computes equilibrium rotational constants from
the current geometry and isotope masses. The Jacobian is evaluated numerically
with respect to the non-redundant GICs, using the GIC B matrix to back-transform
internal-coordinate steps to Cartesian displacements.

Parameters can be frozen with:

```bash
python -m merlino semiexp ... --fixed "GIC001,angle"
```

Each token is matched as a case-insensitive substring of the generated GIC
labels. Fixed parameters are reported but excluded from the least-squares
normal equations.

## Output

The output directory contains:

- `semiexp_geometry.xyz`: fitted equilibrium Cartesian geometry.
- `semiexp_parameters.csv`: final non-redundant GIC values, one-sigma errors and
  active/fixed flags.
- `semiexp_residuals.csv`: observed equilibrium constants, calculated constants
  and residuals in MHz.
- `semiexp_manifest.json`: reproducibility manifest with checksums.

The parameter values use the native Merlino GIC units: stretches in Angstrom and
angular coordinates in radians.
