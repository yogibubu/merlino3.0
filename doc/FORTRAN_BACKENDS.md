# Fortran Backends

Merlino4 keeps only active Fortran backends under `fortran/`. Historical
standalone programs remain in the frozen Merlino3.0 tree.

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

The active build intentionally does not compile old Z-matrix, FITPOT/VCI/DVR,
MSR/isotope or rate utilities. Those historical programs are not duplicated in
Merlino4.

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

## Active Harmonic Internal-Coordinate Source

- Source: `fortran/harmonic_internal/gf.f`
- Compile check:

```bash
cd fortran/harmonic_internal
./compile_check
```

`gf.f` is a source backend rather than a standalone executable. It provides
`DNICGF` for building internal-coordinate F/G matrices and `DNICFq` for Wilson
GF harmonic frequencies and internal-coordinate normal modes. It is the
harmonic internal-coordinate layer that the later anharmonic VPT2/VCI workflow
will build on.

## GDV VPT2/VCI Reference Sources

The active GDV reference sources are inventoried, not copied wholesale:

- `gdv_j32p/gdv/l717.F`: VCI/VPT2 decks, including `VCIDrv`, `VCIGen`,
  `VCIInt`, `VCIPT2` and `VCIVar`.
- `gdv_j32p/gdv/dinautil.F`: anharmonic RWF and harmonic/internal-coordinate
  support, including `AnhFIO`, `DNICGF`, `DNICFq` and `VPT2En`.
- `gdv_j32p/gdv/utilnz.F`: large utility source containing `NHDiag`, the
  Davidson/subspace diagonalization source to extract for large VCI spaces.

See `doc/GDV_VPT2_VCI_SOURCE_MAP.md` before extracting standalone Merlino4
Fortran kernels.

The `fortran/` root is intentionally documentation-only. Generated compiler
logs such as `error` files are not source and should not be committed.
