# Local Symmetry Pipeline

This document summarizes the local-symmetry/topology workflow used to build
non-redundant internal coordinates.

## Overview
- Build atomic equivalence classes from Zeff + topological signatures.
- Recognize local geometries (coord. 4/5/6/8) using Kabsch matching + invariants.
- Group primitives by signatures; build symmetry-local bases.
- Add cyclic ring coordinates (breathing/CVB/cyclic torsions), butterfly, and fused-ring reductions.
- Prune per block with SVD or metric (G) if requested.

## Configuration
```
[u]
mode = auto
symmetry_mode = hybrid
prune_mode = svd
zeff_tol = 0.05
geometry_match_tol = 12.0
pattern_report_path = patterns.json
symmetrize_global = false
keep_a1_only = false
symmetry_tol = 1e-3
symmetry_max_n = 10
assign_symmetry_labels = false
symmetry_quasi_tol = 0.0
symmetry_tol_H = 0.0
heavy_only_orient = false
symmetry_center_idx = -1
ignore_isotopes = false
symmetry_max_dev_strict = 0.0
symmetry_tol_rel = 0.0
symmetry_auto_max_n = false
symmetry_inertia_tol = 1e-3
symmetry_max_radius = 0.0
symmetry_enforce_radial = true
symmetry_profile = false
symmetry_group_limit = c
symmetry_confidence = false
```

Example (advanced global symmetry):
```
[u]
mode = auto
symmetrize_global = true
symmetry_tol_rel = 0.02
symmetry_auto_max_n = true
symmetry_inertia_tol = 1e-3
symmetry_max_radius = 5.0
symmetry_group_limit = c
symmetry_profile = true
symmetry_confidence = true
```

## Cache
Set `MERLINO_FIT_CACHE_DIR` to enable disk caching for `eval_primitives` and `b_matrix`.
This is helpful for large datasets and repeated fits.

## Benchmark
```
python scripts/bench_build_u.py --xyz in.xyz --repeat 5
```
