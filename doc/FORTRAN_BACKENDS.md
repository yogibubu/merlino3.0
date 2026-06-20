# Fortran Backends

Merlino keeps the Fortran code under `fortran/`, separated by backend or
historical module.

## Active PROVA Backend

- Source: `fortran/prova/`
- Main executable in source tree: `fortran/prova/prova`
- Runtime executable used by launchers: `bin/prova.x`
- Build command:

```bash
cd fortran/prova
./compile_MAC
```

`compile_MAC` uses `gfortran -std=legacy -w` because the source is intentionally
fixed-form legacy Fortran. This suppresses compiler deprecation noise such as
shared DO labels, arithmetic IF, and Hollerith constants while preserving
failure on actual compile/link errors.

Build logs are written to `fortran/prova/build/prova_build.log`, which is
ignored by git.

## Other Fortran Areas

- `fortran/gnic/`: standalone/non-active GIC development code.
- `fortran/qcent/`: quadrupole/centering utility backend.
- `fortran/quadrupolari/`: quadrupole conversion examples and utilities.
- `fortran/symmetry/`: standalone symmetry experiments/reference code.
- `fortran/volt/`: reference vibrational input/output material.

Generated compiler logs such as `error` files are not source and should not be
committed.

