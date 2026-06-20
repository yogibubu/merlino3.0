# Harmonic Internal-Coordinate Backend

This directory contains the active Fortran source for harmonic vibrational
analysis in non-redundant internal coordinates.

- `gf.f`: Wilson GF routines in internal coordinates.
  - `DNICGF` builds internal-coordinate F and G matrices from Cartesian
    gradient, Cartesian force constants, B matrix and B-matrix derivatives.
  - `DNICFq` performs the harmonic GF frequency analysis and returns
    frequencies and normal modes in internal coordinates.

`gf.f` is a library source, not a standalone executable. It is compiled as an
object and will be linked later with the shared DiNa/GIC numerical utilities
needed by Merlino4 harmonic and anharmonic workflows.

Use:

```bash
cd fortran/harmonic_internal
./compile_check
```

The check compiles `gf.f` to `build/gf.o` with fixed-form legacy Fortran flags.
The object may still have unresolved external symbols until it is linked into a
complete backend.
