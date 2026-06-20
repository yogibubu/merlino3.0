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
        resolve_source_backend("vpt2_vci", root=root)
        == root / "fortran" / "vpt2_vci" / "vci_core.f"
    )


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


def test_vpt2_vci_fortran_controlled_basis_runtime(tmp_path):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran is not available")
    root = repo_root(Path(__file__))
    driver = tmp_path / "test_basis.f"
    driver.write_text(
        """      Program TestBasis
      Integer Basis(100,3),QMin(3),QMax(3),CMin(4),CMax(4)
      Integer NState,Info,I,Tot,NExc,Bad
      QMin(1)=0
      QMin(2)=0
      QMin(3)=0
      QMax(1)=3
      QMax(2)=2
      QMax(3)=1
      CMin(1)=1
      CMax(1)=2
      CMin(2)=2
      CMax(2)=3
      CMin(3)=3
      CMax(3)=3
      CMin(4)=-1
      CMax(4)=-1
      Call M4VCIBasisCtl(3,4,100,QMin,QMax,CMin,CMax,
     $                   NState,Basis,Info)
      If(Info.ne.0) Stop 10
      If(NState.ne.14) Stop 11
      Bad=0
      Do 30 I=1,NState
         If(Basis(I,2).gt.2) Bad=1
         If(Basis(I,3).gt.1) Bad=1
         Tot=Basis(I,1)+Basis(I,2)+Basis(I,3)
         NExc=0
         If(Basis(I,1).gt.0) NExc=NExc+1
         If(Basis(I,2).gt.0) NExc=NExc+1
         If(Basis(I,3).gt.0) NExc=NExc+1
         If(NExc.eq.1 .and. Tot.gt.2) Bad=1
         If(NExc.eq.2 .and. Tot.gt.3) Bad=1
         If(NExc.eq.3 .and. Tot.ne.3) Bad=1
30    Continue
      If(Bad.ne.0) Stop 12
      End
""",
        encoding="utf-8",
    )
    exe = tmp_path / "test_basis"
    subprocess.run(
        [
            "gfortran",
            "-std=legacy",
            "-ffixed-form",
            str(driver),
            str(root / "fortran" / "vpt2_vci" / "vci_core.f"),
            "-o",
            str(exe),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


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
