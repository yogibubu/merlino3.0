# Merlino Fit Prototype

Independent prototype for robust PES fitting (energies + gradients) with
internal-coordinate basis functions and topology integration.

## Structure
- `survibfit/` core package
- `topology/` copied from Merlino (used for bonding/perception)

## Quick start (prototype)

```
python -m survibfit.cli fit --config example_fit.ini --data data.npz --out terms.txt --with-g
```

## Geometry edit (bond length)
```
python -m survibfit.modify_geom --xyz in.xyz --out out.xyz --bond 1 2 1.10 --xyz-units ang --out-units ang
```

## Geometry edit (BDPCS3 automatic rule)
```
python -m survibfit.modify_geom
```
The CLI will ask for the XYZ path, output path, units, and whether to apply BDPCS3.

Use `--bdpcs3-version updated|legacy` to select the formulation when BDPCS3 is
enabled (default: `updated`). `updated` is the topology/synthon-consistent
BDPCS3 formulation already implemented in Merlino. The back-transform uses
physical metric weights by default: stretches are kept stiff, valence and
out-of-plane coordinates are intermediate, and torsions are softer. Directional
non-covalent H...Y targets are also included for X-H...Y contacts when
X-H-Y >= 150 deg; the correction is damped by an error function in the H...Y
distance so it vanishes for long contacts. Current H-bond values are
O-H...O = -0.055 Angstrom and N-H...N = -0.055 Angstrom; mixed N/O cases use
their average. All H-bond and metric parameters are read from
`merlino_core/parameters/bdpcs3_hbond.toml`; the topology policy is documented
in `doc/HBOND_TOPOLOGY_POLICY.md`. `legacy` is retained only for older runs.

Example (updated BDPCS3):
```
python -m survibfit.modify_geom --bdpcs3-version updated
```

Example (non-interactive updated BDPCS3):
```
python -m survibfit.modify_geom --xyz in.xyz --out out.xyz --rule bdpcs3 --bdpcs3-version updated
```

Add `--report` to print bond changes and rotational constants (GHz).
Use `--isotopes` to select isotopes, e.g. `--isotopes 1:2,3:13`. With `--report` the tool prints the isotopes actually used.

`data.npz` is expected to contain:
- `coords`: shape (ngeom, natoms, 3) in atomic units
- `Z`: shape (natoms,) atomic numbers
- `energy`: shape (ngeom,) in Eh
- `grad`: shape (ngeom, natoms, 3) in Eh/bohr

`example_fit.ini` is provided. `fit.ini` is described in `survibfit/config.py` docstring.

U matrix: you can set `[u] mode = auto` to build a non-redundant internal-coordinate U from the current geometry, or set `u_path` to load a precomputed U.
Additional `[u]` options:
- `symmetry_mode = hybrid|topo|gblock` (default: `hybrid`)
- `prune_mode = svd|g|none` (default: `svd`)
- `zeff_tol = 0.05` (threshold for Zeff-based equivalence classes)
- `geometry_match_tol = 12.0` (deg, local-geometry matching tolerance)
- `pattern_report_path = path/to/patterns.json` (optional report)
- `symmetrize_global = true|false` (default: `false`)
- `keep_a1_only = true|false` (default: `false`)
- `symmetry_tol = 1e-3` (global symmetry tolerance)
- `symmetry_max_n = 10` (max Cn axis tested; include C5 and up to C10)
- `a1_tol = 1e-6` (A1 retention threshold)
- `assign_symmetry_labels = true|false` (emit per-coordinate symmetry labels in report)
- `symmetry_quasi_tol = 0.0` (optional quasi-symmetry tolerance; 0 disables)
- `symmetry_tol_H = 0.0` (optional separate tolerance for H atoms)
- `heavy_only_orient = true|false` (orient by heavy atoms if true)
- `symmetry_center_idx = -1` (optional explicit center atom index)
- `ignore_isotopes = true|false` (treat isotopes as equivalent if true)
- `symmetry_max_dev_strict = 0.0` (optional strict filter on symmetry op deviation)
- `symmetry_tol_rel = 0.0` (relative tolerance factor on matching cutoff)
- `symmetry_auto_max_n = true|false` (reduce max_n for asymmetric inertia)
- `symmetry_inertia_tol = 1e-3` (inertia asymmetry threshold for auto_max_n)
- `symmetry_max_radius = 0.0` (optional scaling radius for matching; 0 uses auto)
- `symmetry_enforce_radial = true|false` (enable/disable radial pre-filter)
- `symmetry_profile = true|false` (include performance stats in report)
- `symmetry_group_limit = c|d|poly` (optional subgroup filter for operations)
- `symmetry_confidence = true|false` (add a confidence score to symmetry report)

## Vibrations (Gaussian log)
```
python -m survibfit.cli vib --log calc.log --fchk calc.fchk --out vib
```
Example scale file: `scale_example.json`

## ReadGIC generator (Gaussian)
```
python -m survibfit.cli gic --xyz in.xyz --out gic.txt
```
By default this command is canonical: it runs GICForge and writes the exact
ReadGIC lines produced by the Fortran backend.  Use `--workdir DIR` to keep the
GICForge working files.  The legacy pure-Python local builder is available only
for diagnostics through `--python-local`; it requires
`MERLINO_ALLOW_PYTHON_LOCAL_GIC=1` and production workflows should not use it
when Python/Fortran identity is required.
The canonical GICForge report lists candidate families before pruning, including
out-of-plane coordinates; only the post-pruning coordinates in `gauin`/`gauin.symm` are part of the
production ReadGIC definition. With `GICSYM`, the promoted `gauin` file contains
symmetry-adapted coordinates and `gicsym` records their irreps.

## Gaussian puckering scan input
Generate a multi-link Gaussian input that defines ring puckering GICs directly
from the current geometry. Rings with at least 4 atoms are supported; tests
cover the 4-, 5-, and 7-membered cases. The writer
emits inactive endocyclic torsions and inactive `RPck...` components; the
effective Gaussian GICs are active `QPck...=SQRT(...)` and
`PhiP...=ATAN2(...)` coordinates. Seven-membered rings generate two
`QPck/PhiP` pairs.

```
python -m survibfit.cli pucker-gaussian \
  --xyz in.xyz \
  --ring 1,2,3,4,5 \
  --phi-start 0 --phi-end 360 --phi-step 10 \
  --gjf-out puckering_scan.gjf \
  --manifest-out puckering_scan.json
```

Use `--constraint-mode functional-targets` for direct phase stepping or
`--constraint-mode scan-to-zero` for the linear fallback on paired puckering
modes.

Options such as `--min-coeff`, `--no-normalize`, `--symmetry-mode` and
`--prune-mode` apply only with `--python-local`.  Canonical GICForge output uses
the Fortran definitions unchanged.
Local symmetry options for `--python-local` (same as `[u]`):
`--symmetry-mode`, `--prune-mode`, `--zeff-tol`, `--geometry-match-tol`, `--pattern-report`,
`--symmetrize-global`, `--keep-a1-only`, `--assign-symmetry-labels`, `--symmetry-quasi-tol`.
`--symmetry-tol-h`, `--heavy-only-orient`.
`--symmetry-center-idx`.
`--ignore-isotopes`.
`--symmetry-max-dev-strict`.
`--symmetry-tol-rel`, `--symmetry-auto-max-n`, `--symmetry-inertia-tol`.
`--symmetry-max-radius`, `--symmetry-no-radial-filter`, `--symmetry-profile`,
`--symmetry-group-limit`, `--symmetry-confidence`.

Example (global symmetrization, keep only A1):
```
python -m survibfit.cli gic --xyz in.xyz --out gic.txt --symmetrize-global --keep-a1-only --assign-symmetry-labels
```

Example (strict symmetry + ignore isotopes):
```
python -m survibfit.cli gic --xyz in.xyz --out gic.txt --symmetrize-global --symmetry-max-dev-strict 1e-4 --ignore-isotopes
```

Example (advanced global options):
```
python -m survibfit.cli gic --xyz in.xyz --out gic.txt --symmetrize-global \
  --symmetry-tol-rel 0.02 --symmetry-auto-max-n --symmetry-inertia-tol 1e-3 \
  --symmetry-max-radius 5.0 --symmetry-group-limit c --symmetry-profile \
  --symmetry-confidence
```

## Molecule similarity (synthons + Gaussian model + ring comparison)
Compare two molecules by fitting a Gaussian model in synthon feature space
(charge, covalency, delocalization, strain, Zeff) and evaluating the
Bhattacharyya distance. The final score can include an additional
ring-aware term (ring size/planarity/position/fusion descriptors):

```
python -m survibfit.synthon_similarity --xyz-a mol1.xyz --xyz-b mol2.xyz
```

Optional outputs and controls:
```
python -m survibfit.synthon_similarity \
  --xyz-a mol1.xyz --xyz-b mol2.xyz \
  --covariance-mode full --regularization 5e-2 --ring-weight 0.25 \
  --json-out similarity.json
```
Use `--no-standardize` to disable global feature standardization.
Use `--no-ring-comparison` to disable the ring-aware term.
Pair mode also prints the top feature contributors (Similarity Explain).

Library mode (one query vs many molecules):
```
python -m survibfit.synthon_similarity \
  --query-xyz query.xyz \
  --library-dir ./library_xyz \
  --library-glob "*.xyz" \
  --top-k 10 \
  --json-out ranking.json
```

Batch mode (many queries vs library) with CSV:
```
python -m survibfit.synthon_similarity \
  --query-dir ./queries_xyz \
  --library-dir ./library_xyz \
  --library-glob "*.xyz" \
  --top-k 5 \
  --json-out batch.json \
  --csv-out batch.csv
```

## Auto-report pipeline
Generate per-query reports including topology report, point-group summary,
and similarity ranking:
```
python -m survibfit.auto_report_pipeline \
  --query-dir ./queries_xyz \
  --library-dir ./library_xyz \
  --out-dir ./auto_reports \
  --top-k 10
```

## Fragment pipeline (MVP)
Fragment-level similarity with GAP detection for library enrichment:
```
python -m survibfit.fragment_pipeline \
  --xyz query.xyz \
  --se-dir /path/to/se_library \
  --pcs2-dir /path/to/pcs2_library \
  --out ./fragment_report \
  --top-k 5 \
  --gap-threshold 0.75
```
PCS2 target mode (use HPCS2 base geometry, PCS2-only matching):
```
python -m survibfit.fragment_pipeline \
  --xyz query.xyz \
  --se-dir /path/to/se_library \
  --pcs2-dir /path/to/pcs2_library \
  --out ./fragment_report \
  --target-pcs2
```
Outputs:
- `fragment_pipeline.json`
- `to_curate.json` (fragments below threshold)
- `fragment_pipeline.md`

PCS2 delta correction from HPCS2 base geometry:
```
python -m survibfit.fragment_delta_correction prepare-hpcs2 \
  --query-xyz query.xyz \
  --fragment-report ./fragment_report/fragment_pipeline.json \
  --hpcs2-dir /path/to/hpcs2_library \
  --out ./fragment_report/delta_bundle_hpcs2

python -m survibfit.fragment_delta_correction apply \
  --query-xyz query.xyz \
  --manifest ./fragment_report/delta_bundle_hpcs2/delta_manifest.json \
  --out-xyz ./fragment_report/query.pcs2.xyz
```

Status:
- the `fragment_delta_correction` module is active code, not an archive helper
- targeted tests live in `tests/test_fragment_delta_correction.py`
- the current repo layout treats `se_library`, `pcs2_library`, and `hpcs2_library`
  as local data libraries rather than source trees

Rovibrational note:
- any local `DeltaVib/alpha` integration in Merlino should be considered a
  downstream compatibility layer
- the primary methodological development for that part now belongs to
  `CeDiTT + alpha_resonances`

## Notes
- Topology perception is performed in Å; coordinates are converted from au internally.
- G and derivatives are computed via finite differences (prototype). B uses analytic gradients for all primitives.
- Gradient/Hessian transforms between Cartesian and internal are available in `survibfit/transforms.py`.
- Coordinate transforms (iterative back-transform) are available in `survibfit/transforms.py`.
- For disconnected systems, TRIC-style fragment coordinates are added: translations use geometric centers and rotations use the quaternion exp-map (Wang & Song, JCP 144, 214108).
- Fortran reference routines are used only for regression testing (optional).
- Legendre basis supports associated orders via `m_i` parameters in config.
- Local symmetry workflow: see `docs/symmetry_local.md`.
