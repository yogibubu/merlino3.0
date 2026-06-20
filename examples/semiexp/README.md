# Semiexperimental Geometry Examples

These examples are small end-to-end inputs for the Merlino Cartesian/GIC
semiexperimental equilibrium-geometry workflow. Each case contains a parent
XYZ file and a TOML isotopologue table accepted by both the CLI and GUI.

Run one case with:

```bash
python -m merlino semiexp \
  --xyz examples/semiexp/water/parent.xyz \
  --observations examples/semiexp/water/isotopologues.toml \
  --outdir working/examples/water_semiexp
```

The outputs include the fitted Cartesian geometry, non-redundant GIC
parameters, propagated errors, residuals, Kraitchman diagnostics, an HTML
report and `semiexp_tables.tex` for manuscript tables.

