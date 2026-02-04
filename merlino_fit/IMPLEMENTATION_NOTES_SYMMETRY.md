# Implementation Notes: Local Symmetry + Topological GIC

## Key Suggestions

1) **Zeff clustering stability**
   - Use `zeff_tol` (default `0.05`) but also quantize Zeff to bins to avoid noise flips.
   - Example: `zbin = round(Zeff / zeff_tol)` and compare bins + |Zeff_i - Zeff_j| < zeff_tol.

2) **Geometry matching confidence**
   - If two local geometries have similar scores (e.g., TBP vs square‑pyramidal),
     mark as ambiguous and fall back to the general basis.
   - Always report ambiguity in the pattern report.

3) **Ring canonicalization**
   - Reuse canonical ring ordering from `RingSet` for all cyclic coordinates.
   - This keeps cyclic bases continuous and deterministic.

4) **Block‑level pruning (SVD)**
   - Prefer SVD per block (angles, cyclic valence bends, cyclic torsions, out-of-plane) to keep interpretability.
   - Global SVD can be used optionally but may mix unrelated coordinate families.

5) **Always generate a pattern report**
   - Even minimal, it should list:
     - local geometry tag + score
     - equivalent classes with Zeff mean
     - chosen symmetry basis (or fallback)
     - ring clusters and reductions applied
     - dihedrals selected per bond

6) **Robust fallback**
   - If any local symmetry cannot be confidently assigned, fall back to the
     generic deterministic basis (sum + orthonormal differences).

7) **Global symmetry acceleration**
   - Cache candidate operations by `max_n`.
   - Match atoms by symbol groups with radial pre-filters and early exits.
   - Optionally reduce `max_n` from inertia asymmetry (auto mode).

8) **Reporting upgrades**
   - Include ring clusters and condensed-ring reductions.
   - Include dihedral selection decisions per bond.
