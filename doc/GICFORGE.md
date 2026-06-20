# GICForge

GICForge is the active Fortran77 geometry backend formerly called `prova`.

## Scope

GICForge is intentionally narrow:

- read Cartesian molecular input from the Merlino working directory
- build topology, redundant GICs and non-redundant GICs
- canonicalize ring atom numbering with the same Prelog-first convention used
  by Python
- optionally build/write the B matrix used for geometry-parameter changes
- write a readable geometry/GIC report
- write Gaussian input with inactive `RPck....` puckering components and active
  `QPck....`/`PhiP....` coordinates

Everything else belongs to Python: GUI orchestration, RDKit/SMILES, project
management, DVR, Cremer-Pople post-processing, regression comparison, freeze
checks and user-facing workflow logic.

## Build

```bash
cd fortran/gicforge
./compile_MAC
```

The build creates:

- `fortran/gicforge/gicforge`
- `bin/gicforge.x`
- `bin/prova.x` as a compatibility alias

The source remains fixed-form Fortran77/legacy-compatible Fortran.

