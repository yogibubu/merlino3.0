# Merlino4 Benchmarks

Small benchmark inputs used to validate the refactor without relying on large
external outputs.

## H2O GF

```bash
python -m merlino gf --fchk gui/tests/gaussian/h2o.fchk --out reports/h2o_gf.txt --csv-dir reports
```

Expected: three positive GF frequencies and a PED table whose columns sum to
100%.

## Toy VPT2/VCI

```bash
python -m merlino vci --qff benchmarks/toy2mode.qff --max-quanta 2 --roots 4 --csv-dir reports
```

Expected: VPT2/VCI comparison report and `vpt2_vci_comparison.csv`.
