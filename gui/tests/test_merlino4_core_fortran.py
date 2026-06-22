from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from merlino_core import ensure_project_state, repo_root, sha256_file, write_manifest
from merlino_fortran import backend_executable, resolve_backend, resolve_source_backend
from merlino_semiexp.fit import (
    _covariance,
    _dynamic_parameter_scales,
    _least_squares_hessian,
    _svd_trust_region_lm_step,
)


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
    assert (
        resolve_source_backend("semiexp", root=root)
        == root / "fortran" / "semiexp" / "semiexp_core.f"
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
    assert (root / "fortran" / "vpt2_vci" / "build" / "vpt2_core.o").exists()


def test_semiexp_fortran_core_runtime(tmp_path):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran is not available")
    root = repo_root(Path(__file__))
    driver = tmp_path / "test_semiexp.f"
    driver.write_text(
        """      Program TestSemiexp
      Integer Info
      Integer ClassMap(2)
      Integer Group(2),NDown
      Double Precision XYZ(3,2),BRow(3,2),Mass(2),ABC(3),PMom(3)
      Double Precision Jac(2,1),Res(2),W(2),DQ(1),Cov(1,1)
      Double Precision Hess(1,1)
      Double Precision JacC(2,2),DQC(2),CovC(2,2),HessC(2,2)
      Double Precision WOut(2),ScaleUsed
      XYZ(1,1)=0.0D0
      XYZ(2,1)=0.0D0
      XYZ(3,1)=0.0D0
      XYZ(1,2)=1.0D0
      XYZ(2,2)=0.0D0
      XYZ(3,2)=0.0D0
      Call M4SEBondB(2,1,2,XYZ,BRow,Info)
      If(Info.ne.0) Stop 10
      If(DAbs(BRow(1,1)+1.0D0).gt.1.0D-12) Stop 11
      If(DAbs(BRow(1,2)-1.0D0).gt.1.0D-12) Stop 12
      Mass(1)=12.0D0
      Mass(2)=16.0D0
      Call M4SERotConst(2,Mass,XYZ,ABC,PMom,Info)
      If(Info.ne.0) Stop 13
      If(ABC(1).le.0.0D0) Stop 14
      Jac(1,1)=1.0D0
      Jac(2,1)=2.0D0
      Res(1)=1.0D0
      Res(2)=1.0D0
      W(1)=1.0D0
      W(2)=1.0D0
      Call M4SENormalEq(2,1,Jac,Res,W,0.0D0,DQ,Cov,Hess,Info)
      If(Info.ne.0) Stop 15
      If(DAbs(DQ(1)-0.6D0).gt.1.0D-10) Stop 16
      If(DAbs(Hess(1,1)-10.0D0).gt.1.0D-10) Stop 17
      JacC(1,1)=1.0D0
      JacC(1,2)=0.0D0
      JacC(2,1)=0.0D0
      JacC(2,2)=1.0D0
      Res(1)=2.0D0
      Res(2)=4.0D0
      ClassMap(1)=1
      ClassMap(2)=1
      Call M4SEClassNormalEq(2,2,1,ClassMap,JacC,Res,W,0.0D0,
     $                       DQC,CovC,HessC,Info)
      If(Info.ne.0) Stop 18
      If(DAbs(DQC(1)-3.0D0).gt.1.0D-10) Stop 19
      If(DAbs(DQC(2)-3.0D0).gt.1.0D-10) Stop 20
      If(DAbs(HessC(1,1)-4.0D0).gt.1.0D-10) Stop 21
      If(DAbs(HessC(1,2)-4.0D0).gt.1.0D-10) Stop 22
      Group(1)=1
      Group(2)=1
      Res(1)=1.0D0
      Res(2)=9.0D0
      Call M4SERobustGroupWeights(2,1,Group,Res,W,1,2.0D0,
     $                            WOut,ScaleUsed,NDown,Info)
      If(Info.ne.0) Stop 23
      If(DAbs(ScaleUsed-2.0D0).gt.1.0D-12) Stop 24
      If(NDown.ne.1) Stop 25
      If(WOut(1).ge.1.0D0.or.WOut(2).ge.1.0D0) Stop 26
      End
""",
        encoding="utf-8",
    )
    exe = tmp_path / "test_semiexp"
    subprocess.run(
        [
            "gfortran",
            "-std=legacy",
            "-ffixed-form",
            str(driver),
            str(root / "fortran" / "semiexp" / "semiexp_core.f"),
            "-o",
            str(exe),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run([str(exe)], check=True, capture_output=True, text=True)


def test_semiexp_fortran_trust_region_matches_python_kernel(tmp_path):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran is not available")
    root = repo_root(Path(__file__))
    driver = tmp_path / "test_semiexp_trust.f"
    driver.write_text(
        """      Program TestSemiexpTrust
      Integer Info,OnBnd,Rank
      Double Precision J(3,2),Res(3),W(3),DQ(2),Cov(2,2),Hess(2,2)
      Double Precision Shift
      J(1,1)=10.0D0
      J(1,2)=0.0D0
      J(2,1)=0.0D0
      J(2,2)=2.0D0
      J(3,1)=1.0D0
      J(3,2)=-1.0D0
      Res(1)=100.0D0
      Res(2)=4.0D0
      Res(3)=5.0D0
      W(1)=1.0D0
      W(2)=0.25D0
      W(3)=2.0D0
      Call M4SETrustNormalEq(3,2,J,Res,W,1.0D-8,0.5D0,DQ,Cov,
     $     Hess,Shift,OnBnd,Rank,Info)
      Write(*,'(9(ES24.16,1X),3I6)') DQ(1),DQ(2),Shift,
     $     Cov(1,1),Cov(1,2),Cov(2,1),Cov(2,2),Hess(1,1),
     $     Hess(1,2),OnBnd,Rank,Info
      End
""",
        encoding="utf-8",
    )
    exe = tmp_path / "test_semiexp_trust"
    subprocess.run(
        [
            "gfortran",
            "-std=legacy",
            "-ffixed-form",
            str(driver),
            str(root / "fortran" / "semiexp" / "semiexp_core.f"),
            "-o",
            str(exe),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    completed = subprocess.run([str(exe)], check=True, capture_output=True, text=True)
    tokens = completed.stdout.split()
    fields = np.asarray([float(item) for item in tokens[:9]], dtype=float)
    ints = [int(item) for item in tokens[9:12]]

    jac = np.array([[10.0, 0.0], [0.0, 2.0], [1.0, -1.0]], dtype=float)
    residual = np.array([100.0, 4.0, 5.0], dtype=float)
    weights = np.array([1.0, 0.25, 2.0], dtype=float)
    sqrt_w = np.sqrt(weights)
    weighted_jac = jac * sqrt_w[:, None]
    weighted_residual = residual * sqrt_w
    scales = _dynamic_parameter_scales(weighted_jac, np.ones(2, dtype=float))
    trust_step = _svd_trust_region_lm_step(
        weighted_jac * scales[None, :],
        weighted_residual,
        damping=1.0e-8,
        trust_radius=0.5,
    )
    expected_dq = scales * trust_step.step
    expected_cov = _covariance(weighted_jac, weighted_residual)
    expected_hess = _least_squares_hessian(weighted_jac)
    singular = np.linalg.svd(weighted_jac * scales[None, :], compute_uv=False)
    tol = max(weighted_jac.shape) * np.finfo(float).eps * max(float(singular[0]), 1.0) * 100.0

    np.testing.assert_allclose(fields[0:2], expected_dq, rtol=2.0e-9, atol=2.0e-9)
    assert fields[2] == pytest.approx(trust_step.shift, rel=2.0e-8, abs=2.0e-8)
    np.testing.assert_allclose(
        np.array([[fields[3], fields[4]], [fields[5], fields[6]]]),
        expected_cov,
        rtol=1.0e-9,
        atol=1.0e-9,
    )
    assert fields[7] == pytest.approx(expected_hess[0, 0], rel=1.0e-12)
    assert fields[8] == pytest.approx(expected_hess[0, 1], rel=1.0e-12)
    assert ints == [1, int(np.sum(singular > tol)), 0]


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


def test_vpt2_fortran_quartic_runtime(tmp_path):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran is not available")
    root = repo_root(Path(__file__))
    driver = tmp_path / "test_vpt2.f"
    driver.write_text(
        """      Program TestVPT2
      Integer Basis(1,1),C3Idx(3,1),C4Idx(4,1)
      Double Precision Freq(1),C3Val(1),C4Val(1)
      Double Precision E0(1),E1(1),E2(1),ETot(1)
      Basis(1,1)=0
      Freq(1)=100.0D0
      C4Idx(1,1)=1
      C4Idx(2,1)=1
      C4Idx(3,1)=1
      C4Idx(4,1)=1
      C4Val(1)=4.0D0
      Call M4VPT2Basis(1,1,Basis,Freq,0,C3Idx,C3Val,
     $                 1,C4Idx,C4Val,E0,E1,E2,ETot)
      If(Abs(E0(1)-50.0D0).gt.1.0D-8) Stop 20
      If(Abs(E1(1)-3.0D0).gt.1.0D-8) Stop 21
      If(Abs(E2(1)).gt.1.0D-8) Stop 22
      If(Abs(ETot(1)-53.0D0).gt.1.0D-8) Stop 23
      End
""",
        encoding="utf-8",
    )
    exe = tmp_path / "test_vpt2"
    subprocess.run(
        [
            "gfortran",
            "-std=legacy",
            "-ffixed-form",
            str(driver),
            str(root / "fortran" / "vpt2_vci" / "vpt2_core.f"),
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
