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

The phthalic anhydride literature-data case uses a job file because it needs
weighted reBO primitive-coordinate predicates:

```bash
python -m merlino semiexp \
  --job examples/semiexp/phthalic_anhydride/phthalic_anhydride_predicates.mse.toml \
  --outdir working/semiexp/phthalic_anhydride
```

The succinic anhydride case is the heavy-atom-only member of the anhydride
series. The production comparison fixes the local hydrogen frame because the
published isotopologues do not substitute H atoms:

```bash
python -m merlino semiexp \
  --job examples/semiexp/succinic_anhydride/succinic_anhydride_fixed_h.mse.toml \
  --outdir working/semiexp/succinic_anhydride_fixed_h
```

The three anhydrides are also available as a multi-molecule class-correction
test.  This fit keeps each molecule's own topology, symmetry and isotopologue
table, but refines shared short/long C-C, carbonyl/single C-O, and C-H stretch
corrections across the homologous series and includes softly regularized bend
classes:

```bash
python -m merlino_core.cli semiexp-ensemble \
  --job examples/semiexp/anhydrides_ensemble/anhydrides_ensemble.mse-ensemble.toml \
  --outdir working/semiexp/anhydrides_ensemble
```

The output includes a text report, CSV summaries, covariance/correlation
matrices and `ensemble_manifest.json`.  Ensemble reports include an explicit
model status: `accepted` means the model passes the rank, conditioning and
support checks; `review` means the fit is usable but contains flagged class
correlations; `rejected` means the model should not be used for production
corrections without changing classes, priors or atom typing.

The complete no-prior/soft-prior/hard-constraint comparison, prior-strength
scan, leave-one-molecule-out diagnostics, and JPCL LaTeX fragments can be
regenerated with:

```bash
python -m merlino_core.cli semiexp-ensemble-paper \
  --job examples/semiexp/anhydrides_ensemble/anhydrides_ensemble.mse-ensemble.toml \
  --paper-dir doc/papers/ensemble_jpcl \
  --outdir working/semiexp/anhydrides_full_analysis
```

The two glycine conformers can be fitted together directly from the legacy MSR
benchmark inputs.  The `synthon` variant uses continuous effective atomic
numbers (`Zeff`) with a threshold to define transferable atom types across the
two different conformer numberings:

```bash
python -m merlino_core.cli semiexp-ensemble \
  --job examples/semiexp/glycine_ensemble/glycine_conformers_synthon.mse-ensemble.toml \
  --outdir working/semiexp/glycine_ensemble_synthon
```

Scan the continuous atom-typing threshold with:

```bash
python -m merlino_core.cli semiexp-ensemble-synthon-scan \
  --job examples/semiexp/glycine_ensemble/glycine_conformers_synthon.mse-ensemble.toml \
  --outdir working/semiexp/glycine_synthon_threshold_scan_wide \
  --threshold 0.010 --threshold 0.035 --threshold 0.075 \
  --threshold 0.100 --threshold 0.150 --threshold 0.250
```
