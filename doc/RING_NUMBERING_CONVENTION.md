# Ring Numbering Convention

Merlino uses a deterministic cyclic numbering before building ring-puckering
GICs.

## Canonical Rule

For every detected ring:

1. Keep only a cyclic sequence of bonded ring atoms.
2. Rotate the sequence so the lowest input atom index is first.
3. Choose the direction so the second atom is the lower of the two neighbours
   of the first atom.

Example:

```text
3 4 5 1 2  ->  1 2 3 4 5
4 3 2 1 5  ->  1 2 3 4 5
```

This removes arbitrary DFS/RDKit traversal choices. For perfectly symmetric
rings, the absolute phase origin remains conventional; this rule fixes that
origin reproducibly from the input atom numbering.

## Python

`merlino_fit.survibfit.puckering_gaussian` canonicalizes both user-provided
`--ring` sequences and automatically detected rings before writing:

- cyclic torsions `T...`
- inactive puckering components `RPck...`
- active Gaussian coordinates `QPck...` and `PhiP...`

The DVR Cremer-Pople labeling uses the same canonical sequence.

## Fortran

The Fortran `prova` path already canonicalizes cycles in `mkcyc.f`:

- `CanCyc` rotates the cycle to the canonical atom order.
- `SymCyc` is intentionally not applied before the ring-puckering GIC build,
  because a symmetry-based shift would change the phase origin.
- `CyGND` builds ring puckering dihedrals from the canonical cycle.

Therefore Fortran and Python now share the same deterministic ring order. If
`IPrint > 0`, `prova` prints the input and canonical cycle numbering so
mismatches can be diagnosed directly.

## Practical Implication

Changing the atom numbering in the input can change the reported phase origin,
but rerunning the same input or reversing/rotating a detected cycle will not.
`Q` is invariant to this convention; `Phi` is convention-dependent and should
always be interpreted together with the documented ring order.
