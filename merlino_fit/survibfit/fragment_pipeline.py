from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .modify_geom import read_xyz
from .pipeline import _load_topology_elements, build_topology_full


@dataclass
class Fragment:
    source_xyz: str
    source_library: str
    molecule_key: str
    fragment_id: str
    atom_indices: list[int]
    feature_vector: np.ndarray
    ring_size: int


def _xyz_to_au(path: Path):
    atoms, coords_ang, _ = read_xyz(path)
    atomic_number = _load_topology_elements()
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords_au = np.array(coords_ang, dtype=float) / 0.52917721092
    return atoms, Z, coords_au


def _ring_fragment_features(path: Path, source_library: str) -> list[Fragment]:
    atoms, Z, coords_au = _xyz_to_au(path)
    _cg, _dg, ringset, synthons, _arom = build_topology_full(coords_au, Z)
    if ringset is None or len(ringset) == 0:
        return []

    coords = np.array(coords_au, dtype=float)
    mol_center = coords.mean(axis=0)
    mol_scale = float(np.sqrt(((coords - mol_center) ** 2).sum(axis=1).mean()))
    if mol_scale < 1.0e-12:
        mol_scale = 1.0

    out = []
    for idx, ring in enumerate(ringset.rings):
        ring_atoms = sorted(set(int(a) for a in ring.atoms))
        # MVP neighborhood: ring atoms + first-shell neighbors.
        neigh = set(ring_atoms)
        for a in ring_atoms:
            for b in range(len(atoms)):
                if b == a:
                    continue
                # geometric first-shell using BO > 0.3 from topology model
                # via synthons helper bond_order on atom pair.
                bo = synthons.bond_order(a, b)
                if bo > 0.3:
                    neigh.add(int(b))
        frag_atoms = sorted(neigh)

        syn = np.zeros(5, dtype=float)
        for a in frag_atoms:
            syn[0] += float(synthons.charge(a))
            syn[1] += float(synthons.covalency(a))
            syn[2] += float(synthons.delocalization(a))
            syn[3] += float(synthons.strain(a))
            syn[4] += float(synthons.Zeff(a))
        syn /= max(len(frag_atoms), 1)

        xyz_ring = coords[ring_atoms]
        ctr = xyz_ring.mean(axis=0)
        ring_size = float(len(ring_atoms))
        planarity = float(ring.planarity() or 0.0)
        fused_degree = float(len(ring.connected_rings))
        radial_pos = float(np.linalg.norm(ctr - mol_center) / mol_scale)
        ring_radius = float(np.sqrt(((xyz_ring - ctr) ** 2).sum(axis=1).mean()) / mol_scale)
        frag_size = float(len(frag_atoms))

        vec = np.concatenate(
            [
                syn,
                np.array(
                    [ring_size, planarity, fused_degree, radial_pos, ring_radius, frag_size],
                    dtype=float,
                ),
            ]
        )
        out.append(
            Fragment(
                source_xyz=str(path),
                source_library=source_library,
                molecule_key=path.stem.lower(),
                fragment_id=f"{path.stem}:ring_{idx+1}",
                atom_indices=frag_atoms,
                feature_vector=vec,
                ring_size=int(ring_size),
            )
        )
    return out


def _sim(a: np.ndarray, b: np.ndarray) -> float:
    # Stable bounded similarity in (0,1].
    d = float(np.linalg.norm(a - b))
    return float(np.exp(-d))


def _rank_fragment(query_frag: Fragment, library_frags: list[Fragment], top_k: int):
    rows = []
    for cand in library_frags:
        if cand.source_xyz == query_frag.source_xyz and cand.fragment_id == query_frag.fragment_id:
            continue
        s = _sim(query_frag.feature_vector, cand.feature_vector)
        rows.append(
            {
                "candidate_xyz": cand.source_xyz,
                "candidate_library": cand.source_library,
                "candidate_molecule_key": cand.molecule_key,
                "candidate_fragment_id": cand.fragment_id,
                "similarity": s,
                "candidate_ring_size": cand.ring_size,
            }
        )
    # Prefer SE over PCS2 when the same molecule exists in both libraries.
    by_molecule = {}
    for row in rows:
        key = row["candidate_molecule_key"]
        cur = by_molecule.get(key)
        if cur is None:
            by_molecule[key] = row
            continue
        # Same molecule key: keep SE candidate if available, otherwise best score.
        if cur["candidate_library"] != "SE" and row["candidate_library"] == "SE":
            by_molecule[key] = row
            continue
        if cur["candidate_library"] == row["candidate_library"]:
            if row["similarity"] > cur["similarity"]:
                by_molecule[key] = row

    uniq_rows = list(by_molecule.values())
    uniq_rows.sort(
        key=lambda r: (r["similarity"], 1 if r["candidate_library"] == "SE" else 0),
        reverse=True,
    )
    return uniq_rows[: max(1, int(top_k))]


def run_fragment_pipeline(
    xyz: Path,
    se_dir: Path,
    pcs2_dir: Path,
    out_dir: Path,
    *,
    xyz_glob: str = "*.xyz",
    top_k: int = 5,
    gap_threshold: float = 0.75,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    q_frags = _ring_fragment_features(xyz, source_library="QUERY")

    lib_xyz = sorted(se_dir.glob(xyz_glob)) + sorted(pcs2_dir.glob(xyz_glob))
    lib_frags = []
    for p in sorted(se_dir.glob(xyz_glob)):
        try:
            lib_frags.extend(_ring_fragment_features(p, source_library="SE"))
        except Exception:
            continue
    for p in sorted(pcs2_dir.glob(xyz_glob)):
        try:
            lib_frags.extend(_ring_fragment_features(p, source_library="PCS2"))
        except Exception:
            continue

    report_rows = []
    to_curate = []
    for qf in q_frags:
        ranking = _rank_fragment(qf, lib_frags, top_k=top_k)
        best = ranking[0] if ranking else None
        best_sim = float(best["similarity"]) if best is not None else 0.0
        status = "OK" if best_sim >= gap_threshold else "GAP"
        row = {
            "query_fragment_id": qf.fragment_id,
            "query_xyz": qf.source_xyz,
            "query_ring_size": qf.ring_size,
            "best_similarity": best_sim,
            "status": status,
            "ranking": ranking,
        }
        report_rows.append(row)
        if status == "GAP":
            to_curate.append(
                {
                    "query_fragment_id": qf.fragment_id,
                    "query_xyz": qf.source_xyz,
                    "best_similarity": best_sim,
                    "best_candidate": best,
                    "action": "Add a closer high-accuracy fragment to SE/PCS2 library.",
                }
            )

    report = {
        "query_xyz": str(xyz),
        "query_fragments": len(q_frags),
        "library_xyz_count": len(lib_xyz),
        "library_fragments": len(lib_frags),
        "gap_threshold": float(gap_threshold),
        "fragments": report_rows,
        "to_curate": to_curate,
    }

    json_path = out_dir / "fragment_pipeline.json"
    curate_path = out_dir / "to_curate.json"
    md_path = out_dir / "fragment_pipeline.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    curate_path.write_text(json.dumps(to_curate, indent=2) + "\n", encoding="utf-8")

    lines = [
        f"# Fragment Pipeline Report: {xyz.name}",
        "",
        f"- Query fragments: `{len(q_frags)}`",
        f"- Library XYZ count: `{len(lib_xyz)}`",
        f"- Library fragments: `{len(lib_frags)}`",
        f"- Gap threshold: `{gap_threshold:.2f}`",
        "",
        "## Fragment Matches",
        "",
        "| Fragment | Ring size | Best similarity | Status |",
        "|---|---:|---:|---|",
    ]
    for r in report_rows:
        lines.append(
            f"| `{r['query_fragment_id']}` | {r['query_ring_size']} | "
            f"{r['best_similarity']:.6f} | {r['status']} |"
        )
    lines.append("")
    lines.append(f"GAP fragments: `{len(to_curate)}`")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fragment similarity + GAP detection MVP.")
    ap.add_argument("--xyz", required=True, help="Query molecule XYZ.")
    ap.add_argument("--se-dir", required=True, help="SE library directory.")
    ap.add_argument("--pcs2-dir", required=True, help="PCS2 library directory.")
    ap.add_argument("--out", required=True, help="Output directory.")
    ap.add_argument("--xyz-glob", default="*.xyz")
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--gap-threshold", type=float, default=0.75)
    args = ap.parse_args(argv)

    rep = run_fragment_pipeline(
        Path(args.xyz),
        Path(args.se_dir),
        Path(args.pcs2_dir),
        Path(args.out),
        xyz_glob=args.xyz_glob,
        top_k=args.top_k,
        gap_threshold=args.gap_threshold,
    )
    print(f"Fragment report written to: {Path(args.out) / 'fragment_pipeline.json'}")
    print(f"Fragments: {rep['query_fragments']}  GAPs: {len(rep['to_curate'])}")


if __name__ == "__main__":
    main()
