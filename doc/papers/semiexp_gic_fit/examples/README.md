# Semiexperimental MSR Compatibility Examples

`nitrobenzene_cartesian_constraints.msr` is the Cartesian/ModRedundant
equivalent of `nitrobenzene_all_msr.inp`.

The original MSR Z-matrix was converted to Cartesian coordinates with the
Merlino MSR reader. Z-matrix variables marked with `#` are represented as
Gaussian-style ModRedundant freeze records in the `constraints` block. The
isotopologue mass blocks and rotational-data sections are preserved in MSR
compatibility format.
The same `constraints` block can also contain Gaussian-style GIC expression
constraints, for example `QFIX=[GIC001+2*GIC002] Value=0.0` or
`DR(Frozen,Value=0.0)=R[1,3]-R[1,2]`. The Gaussian keyword forms
`Label(Options)=Expression` and `Expression Options` are accepted for hard
constraints, for example `HOH(Frozen)=A(2,1,3)` or
`DR(Frozen,Value=0.0)=R[1,3]-R[1,2]`.

The file can be checked with:

```bash
python -m merlino semiexp \
  --job doc/papers/semiexp_gic_fit/examples/nitrobenzene_cartesian_constraints.msr \
  --outdir working/semiexp/nitrobenzene/run_cartesian_constraints \
  --coordinate-model gic \
  --rotational-components auto
```
