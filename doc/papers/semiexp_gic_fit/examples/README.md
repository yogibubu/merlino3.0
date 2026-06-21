# Semiexperimental MSR Compatibility Examples

`nitrobenzene_cartesian_constraints.msr` is the Cartesian/ModRedundant
equivalent of `nitrobenzene_all_msr.inp`.

The original MSR Z-matrix was converted to Cartesian coordinates with the
Merlino MSR reader. Z-matrix variables marked with `#` are represented as
Gaussian-style ModRedundant freeze records in the `constraints` block. The
isotopologue mass blocks and rotational-data sections are preserved in MSR
compatibility format.

The file can be checked with:

```bash
python -m merlino semiexp \
  --job doc/papers/semiexp_gic_fit/examples/nitrobenzene_cartesian_constraints.msr \
  --outdir working/semiexp/nitrobenzene/run_cartesian_constraints \
  --coordinate-model gic \
  --rotational-components auto
```
