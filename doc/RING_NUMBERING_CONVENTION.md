# Ring Numbering Convention

Merlino uses a deterministic cyclic numbering before building ring-puckering
GICs.

## Prelog-First Canonical Rule

For every detected ring:

1. Keep only a cyclic sequence of bonded ring atoms.
2. Choose the starting atom from local Prelog/CIP priority. Atomic number is
   primary; local connectivity and exocyclic substituents refine the ordering
   when the molecular graph is available.
3. Choose the direction by comparing the two cyclic paths from that atom with
   the same priority rule.
4. Use the input atom index only as the final deterministic tie-break for truly
   equivalent atoms.

Example:

```text
C C C O C  ->  O C C C C
3 4 5 1 2  ->  Prelog start, then deterministic cyclic direction
```

This removes arbitrary DFS/RDKit traversal choices while keeping the chemically
meaningful Prelog origin. For perfectly symmetric rings, the absolute phase
origin remains conventional; the final input-index tie-break fixes that origin
reproducibly.

## Python

`merlino_fit.survibfit.puckering_gaussian` canonicalizes both user-provided
`--ring` sequences and automatically detected rings before writing:

- cyclic torsions `T...`
- inactive puckering components `RPck...`
- active Gaussian coordinates `QPck...` and `PhiP...`

The DVR Cremer-Pople labeling uses the same canonical sequence.

## Fortran

The Fortran `GICForge` path already canonicalizes cycles in `mkcyc.f`:

- `CanCyc` rebuilds a bonded cycle and then applies the same Prelog-first
  rotation/direction rule used by Python.
- `SymCyc` is intentionally not applied before the ring-puckering GIC build,
  because a symmetry-based shift would change the phase origin.
- `CyGND` builds ring puckering dihedrals from the canonical cycle.

Therefore Fortran and Python now share the same deterministic ring order. If
`IPrint > 0`, GICForge prints the input and canonical cycle numbering so
mismatches can be diagnosed directly.

## Practical Implication

Changing the atom numbering in the input can change the reported phase origin,
but rerunning the same input or reversing/rotating a detected cycle will not.
`Q` is invariant to this convention; `Phi` is convention-dependent and should
always be interpreted together with the documented ring order.
