# Fortran Backends

Merlino keeps the Fortran code under `fortran/`, separated by backend or
historical module.

## Active GICForge Backend

- Source: `fortran/gicforge/`
- Runtime executable used by launchers: `bin/gicforge.x`
- Compatibility runtime alias: `bin/prova.x`
- Build command:

```bash
cd fortran/gicforge
./compile_MAC
```

`compile_MAC` uses `gfortran -std=legacy -w` because the source is intentionally
fixed-form legacy Fortran. This suppresses compiler deprecation noise such as
shared DO labels, arithmetic IF, and Hollerith constants while preserving
failure on actual compile/link errors.

The source-tree build output is `fortran/gicforge/build/gicforge`; runtime
copies are written to `bin/gicforge.x` and `bin/prova.x`. Build logs are written
to `fortran/gicforge/build/gicforge_build.log`. The whole `build/` directory is
ignored by git.

GICForge receives only Cartesian XYZ input from Merlino, builds redundant and
non-redundant GICs, optionally writes the B matrix, emits a readable report, and
creates Gaussian input. GUI orchestration, RDKit/SMILES, any conversion to
Cartesian coordinates, DVR and post-processing remain Python responsibilities.

The active build intentionally does not compile the old Z-matrix,
FITPOT/VCI/DVR, MSR/isotope or rate utilities. Historical vibrational material
is kept under `fortran/legacy/`; old non-Cartesian input experiments are kept
outside the active Fortran tree under `archives/legacy_fortran/`.

## Active DVR Backend

- Source: `fortran/dvr/`
- Runtime executable used by launchers: `bin/path_dvr.x`
- Build command:

```bash
cd fortran/dvr
./compile_MAC
```

The Fortran77 DVR backend reads only `dvrin`. Python remains responsible for
reading Gaussian output, constructing path/grid CSV files and launching the
bridge script. The backend currently supports one-dimensional grid DVR,
one-dimensional distributed Gaussian basis DVR and a two-dimensional product
grid DVR.

Diagonalization is performed by `DVRHQRII` in `fortran/dvr/dvr_hqrii.f`, a
renamed local copy of GICForge `HQRII1`. This avoids Jacobi diagonalization for
large Hamiltonians.

## Other Fortran Areas

- `fortran/gnic/`: standalone/non-active GIC development code.
- `fortran/legacy/`: historical vibrational and reader utilities no longer
  compiled into active backends.
- `fortran/qcent/`: quadrupole/centering utility backend.
- `fortran/quadrupolari/`: quadrupole conversion examples and utilities.
- `fortran/symmetry/`: standalone symmetry experiments/reference code.
- `fortran/volt/`: reference vibrational input/output material.

The `fortran/` root is intentionally documentation-only. Generated compiler
logs such as `error` files are not source and should not be committed.
