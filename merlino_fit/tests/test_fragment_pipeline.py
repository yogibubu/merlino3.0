from pathlib import Path
import shutil

from survibfit.fragment_pipeline import run_fragment_pipeline


def test_fragment_pipeline_outputs_reports(tmp_path):
    root = Path(__file__).resolve().parent / "data"
    out = run_fragment_pipeline(
        root / "naphthalene_c10.xyz",
        root,
        root,
        tmp_path,
        top_k=3,
        gap_threshold=0.8,
    )
    assert "fragments" in out
    assert "to_curate" in out
    assert (tmp_path / "fragment_pipeline.json").exists()
    assert (tmp_path / "to_curate.json").exists()
    assert (tmp_path / "fragment_pipeline.md").exists()


def test_fragment_pipeline_prefers_se_for_same_molecule(tmp_path):
    root = Path(__file__).resolve().parent / "data"
    query = root / "naphthalene_c10.xyz"
    se_dir = tmp_path / "se"
    pcs2_dir = tmp_path / "pcs2"
    out_dir = tmp_path / "out"
    se_dir.mkdir()
    pcs2_dir.mkdir()

    # Same molecule key in both libraries -> preference must go to SE.
    shutil.copy2(query, se_dir / "samemol.xyz")
    shutil.copy2(query, pcs2_dir / "samemol.xyz")

    out = run_fragment_pipeline(
        query,
        se_dir,
        pcs2_dir,
        out_dir,
        top_k=1,
        gap_threshold=0.5,
    )
    assert out["fragments"]
    for row in out["fragments"]:
        if row["ranking"]:
            assert row["ranking"][0]["candidate_library"] == "SE"
