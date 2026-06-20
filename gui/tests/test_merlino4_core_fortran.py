from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from merlino_core import ensure_project_state, repo_root, sha256_file, write_manifest
from merlino_fortran import backend_executable, resolve_backend, resolve_source_backend


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


def test_fortran_source_backend_resolution():
    root = repo_root(Path(__file__))
    assert (
        resolve_source_backend("harmonic_internal", root=root)
        == root / "fortran" / "harmonic_internal" / "gf.f"
    )
    assert (
        resolve_source_backend("vpt2_vci", root=root)
        == root / "fortran" / "vpt2_vci" / "vci_core.f"
    )


def test_harmonic_internal_source_compiles_to_object():
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran is not available")
    root = repo_root(Path(__file__))
    subprocess.run(
        [str(root / "fortran" / "harmonic_internal" / "compile_check")],
        cwd=root / "fortran" / "harmonic_internal",
        check=True,
        capture_output=True,
        text=True,
    )
    assert (root / "fortran" / "harmonic_internal" / "build" / "gf.o").exists()


def test_vpt2_vci_source_compiles_to_objects():
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran is not available")
    root = repo_root(Path(__file__))
    subprocess.run(
        [str(root / "fortran" / "vpt2_vci" / "compile_check")],
        cwd=root / "fortran" / "vpt2_vci",
        check=True,
        capture_output=True,
        text=True,
    )
    assert (root / "fortran" / "vpt2_vci" / "build" / "gf_core.o").exists()
    assert (root / "fortran" / "vpt2_vci" / "build" / "vci_core.o").exists()
    assert (root / "fortran" / "vpt2_vci" / "build" / "davidson_core.o").exists()


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
