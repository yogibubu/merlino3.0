from pathlib import Path
import json

import numpy as np

from survibfit.fragment_delta_correction import (
    prepare_delta_workflow,
    prepare_low_level_gaussian_inputs,
    prepare_high_level_curation_jobs,
    apply_delta_correction,
)
from survibfit.modify_geom import write_xyz, read_xyz


def test_prepare_delta_workflow_builds_manifest(tmp_path):
    query = tmp_path / "query.xyz"
    cand = tmp_path / "cand.xyz"
    atoms = ["C", "N", "O", "H"]
    q = np.array(
        [[0.0, 0.0, 0.0], [1.3, 0.1, 0.0], [0.2, 1.2, 0.0], [2.0, 2.0, 0.0]],
        dtype=float,
    )
    c = q + np.array([0.05, -0.02, 0.01])
    write_xyz(query, atoms, q, comment="q")
    write_xyz(cand, atoms, c, comment="c")

    report = {
        "fragments": [
            {
                "query_fragment_id": "q:ring_1",
                "query_atom_indices": [0, 1, 2],
                "ranking": [
                    {
                        "candidate_xyz": str(cand),
                        "candidate_atom_indices": [0, 1, 2],
                        "candidate_fragment_label": "cand:ring_1 (cycle_5)",
                        "similarity": 0.8,
                    }
                ],
            }
        ]
    }
    report_path = tmp_path / "fragment_pipeline.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")

    out = tmp_path / "delta_bundle"
    manifest = prepare_delta_workflow(query, report_path, out)
    assert len(manifest["entries"]) == 1
    e0 = manifest["entries"][0]
    assert e0["high_level_result_xyz"].endswith("candidate_fragment.xyz")
    assert e0["status"] == "pending_low_level_only"
    assert "suggested_new_fragment" in e0
    assert Path(e0["proposed_new_library_fragment_xyz"]).exists()
    assert (out / "delta_manifest.json").exists()
    frag_dir = next(out.glob("frag_*"))
    assert (frag_dir / "query_fragment.xyz").exists()
    assert (frag_dir / "candidate_fragment.xyz").exists()


def test_apply_delta_correction_transfers_local_deformation(tmp_path):
    query = tmp_path / "query.xyz"
    atoms = ["C", "N", "O", "H"]
    q = np.array(
        [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.1, 1.0, 0.0], [3.0, 3.0, 0.0]],
        dtype=float,
    )
    write_xyz(query, atoms, q, comment="q")

    frag_dir = tmp_path / "frag"
    frag_dir.mkdir()
    low = frag_dir / "low_level_result.xyz"
    high = frag_dir / "high_level_result.xyz"
    low_xyz = q[:3].copy()
    high_xyz = low_xyz.copy()
    high_xyz[0, 0] += 0.08
    high_xyz[1, 1] -= 0.06
    high_xyz[2, 0] -= 0.04
    write_xyz(low, atoms[:3], low_xyz, comment="low")
    write_xyz(high, atoms[:3], high_xyz, comment="high")

    manifest = {
        "entries": [
            {
                "query_fragment_id": "q:ring_1",
                "query_atom_indices": [0, 1, 2],
                "low_level_result_xyz": str(low),
                "high_level_result_xyz": str(high),
                "weight_prior": 1.0,
            }
        ]
    }
    manifest_path = tmp_path / "delta_manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    out_xyz = tmp_path / "corrected.xyz"
    meta = apply_delta_correction(query, manifest_path, out_xyz)
    assert meta["entries_used"] == 1
    out_atoms, out_coords, _ = read_xyz(out_xyz)
    assert out_atoms == atoms
    # Fragment atoms should move non-trivially.
    moved = np.linalg.norm(out_coords[:3] - q[:3], axis=1)
    assert float(np.max(moved)) > 1.0e-3
    # Non-fragment atom should remain unchanged.
    assert np.allclose(out_coords[3], q[3], atol=1.0e-8)


def test_prepare_low_level_gaussian_inputs_writes_gjf(tmp_path):
    query = tmp_path / "query.xyz"
    cand = tmp_path / "cand.xyz"
    atoms = ["C", "N", "O"]
    q = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.1, 1.0, 0.0]], dtype=float)
    c = q + 0.03
    write_xyz(query, atoms, q, comment="q")
    write_xyz(cand, atoms, c, comment="c")

    report = {
        "fragments": [
            {
                "query_fragment_id": "q:ring_1",
                "query_atom_indices": [0, 1, 2],
                "ranking": [
                    {
                        "candidate_xyz": str(cand),
                        "candidate_atom_indices": [0, 1, 2],
                        "candidate_fragment_label": "cand:ring_1 (cycle_5)",
                        "similarity": 0.8,
                    }
                ],
            }
        ]
    }
    report_path = tmp_path / "fragment_pipeline.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    bundle = tmp_path / "bundle"
    prepare_delta_workflow(query, report_path, bundle)
    manifest_path = bundle / "delta_manifest.json"

    out = prepare_low_level_gaussian_inputs(
        manifest_path,
        route="#p hf/3-21g opt",
        charge=0,
        multiplicity=1,
        nproc=4,
        mem="2GB",
    )
    assert out["generated_inputs"] == 1
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    e0 = manifest["entries"][0]
    gjf = Path(e0["low_level_gaussian_input"])
    assert gjf.exists()
    txt = gjf.read_text(encoding="utf-8")
    assert "%nprocshared=4" in txt
    assert "%mem=2GB" in txt
    assert "#p hf/3-21g opt" in txt


def test_prepare_high_level_curation_jobs_writes_queue(tmp_path):
    query = tmp_path / "query.xyz"
    cand = tmp_path / "cand.xyz"
    atoms = ["C", "N", "O"]
    q = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.1, 1.0, 0.0]], dtype=float)
    c = q + 0.03
    write_xyz(query, atoms, q, comment="q")
    write_xyz(cand, atoms, c, comment="c")

    report = {
        "fragments": [
            {
                "query_fragment_id": "q:ring_1",
                "query_atom_indices": [0, 1, 2],
                "is_low_score": True,
                "suggested_new_fragment": {
                    "closest_library_reference": {
                        "candidate_fragment_label": "cand:ring_1 (ring_3)"
                    }
                },
                "ranking": [
                    {
                        "candidate_xyz": str(cand),
                        "candidate_atom_indices": [0, 1, 2],
                        "candidate_fragment_label": "cand:ring_1 (cycle_5)",
                        "similarity": 0.3,
                    }
                ],
            }
        ]
    }
    report_path = tmp_path / "fragment_pipeline.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    bundle = tmp_path / "bundle"
    prepare_delta_workflow(query, report_path, bundle)
    out = prepare_high_level_curation_jobs(
        bundle / "delta_manifest.json",
        tmp_path / "hl_jobs",
        route="#p wb97xd/def2tzvp opt",
        charge=0,
        multiplicity=1,
        nproc=4,
        mem="8GB",
    )
    assert len(out["jobs"]) == 1
    gjf = Path(out["jobs"][0]["high_level_input_gjf"])
    assert gjf.exists()
    txt = gjf.read_text(encoding="utf-8")
    assert "#p wb97xd/def2tzvp opt" in txt
