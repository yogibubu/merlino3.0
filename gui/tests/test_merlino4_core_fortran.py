from __future__ import annotations

import json
from pathlib import Path

from merlino_core import ensure_project_state, repo_root, sha256_file, write_manifest
from merlino_fortran import backend_executable, resolve_backend


def test_repo_root_finds_merlino4_root():
    root = repo_root(Path(__file__))
    assert (root / "README.md").exists()
    assert (root / "doc" / "MERLINO4_REFACTOR_PLAN.md").exists()


def test_manifest_helpers_are_deterministic(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("merlino4\n", encoding="utf-8")
    manifest = write_manifest(
        tmp_path / "manifest.json",
        {"sha256": sha256_file(source), "name": "test"},
    )
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert data == {
        "name": "test",
        "sha256": "5bfcef976def954e5384b7d8e0034783edb461fc5f07b63546a8e026de694e58",
    }


def test_fortran_backend_resolution_uses_repo_bin():
    root = repo_root(Path(__file__))
    assert backend_executable("gicforge") == "gicforge.x"
    assert resolve_backend("gicforge", root=root) == root / "bin" / "gicforge.x"
    assert resolve_backend("dvr", root=root) == root / "bin" / "path_dvr.x"


def test_project_state_is_created_and_reloaded(tmp_path):
    state = ensure_project_state(tmp_path)
    assert state.workdir == tmp_path
    assert state.active_workflow == "molecule"
    assert (tmp_path / ".merlino" / "project.json").exists()

    state.active_workflow = "dvr"
    state.metadata["note"] = "roundtrip"
    state.save()

    loaded = ensure_project_state(tmp_path)
    assert loaded.active_workflow == "dvr"
    assert loaded.metadata == {"note": "roundtrip"}
