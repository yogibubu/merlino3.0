# GICForge

GICForge is the active Fortran77 geometry backend formerly called `prova`.

## Scope

GICForge is intentionally narrow:

- read Cartesian molecular input from the Merlino working directory
- build topology, redundant GICs and non-redundant GICs
- remove residual linear dependencies by coordinate-type blocks, so stretches,
  bends, linear bends, torsions and out-of-plane coordinates are never mixed
  during the final rank pruning
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

## Residual Redundancy Pruning

After the ordinary GIC construction, GICForge calls `PruneGICBlocks` before
writing Gaussian input. The routine builds the B matrix for the current GIC
candidates and applies a modified Gram-Schmidt rank test separately to each
coordinate family:

- stretch
- bend
- linear bend
- torsion, including ring puckering and butterfly coordinates
- out-of-plane

Dependent rows are removed from their own family and the remaining arrays are
compacted in place. This is intentionally analogous to the Python
`prune_mode=svd/g` policy, but preserves strict Fortran77 implementation and
does not combine heterogeneous coordinate types.
