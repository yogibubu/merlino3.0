from pathlib import Path

import numpy as np
import pytest

from gui.bdpcs3_workflow import compute_bdpcs3, write_bdpcs3_outputs


def _write_min_xyzin(path: Path):
    text = """3
H2O test
O 0.000000 0.000000 0.000000
H 0.757160 0.586260 0.000000
H -0.757160 0.586260 0.000000
#BASIC
CHARGE              0
SPIN_MULTIPLICITY   1
POINT_GROUP         C1
REPRESENTATION      Ir
T_K =               298.15
P_ATM =             1.000000
"""
    path.write_text(text)


def _write_hbond_xyzin(path: Path):
    text = """3
OH...O test
O 0.000000 0.000000 0.000000
H 0.960000 0.000000 0.000000
O 2.772039 0.319493 0.000000
#BASIC
CHARGE              0
SPIN_MULTIPLICITY   1
POINT_GROUP         C1
REPRESENTATION      Ir
T_K =               298.15
P_ATM =             1.000000
"""
    path.write_text(text)


def test_bdpcs3_workflow_outputs(tmp_path):
    xyzin = tmp_path / "xyzin"
    _write_min_xyzin(xyzin)

    result = compute_bdpcs3(xyzin, "legacy")
    assert len(result.symbols) == 3
    assert result.coords_dpcs3_ang.shape == (3, 3)
    assert result.coords_bdpcs3_ang.shape == (3, 3)
    assert np.isfinite(result.B0_mhz).all()
    assert np.isfinite(result.B1_mhz).all()

    report_path, xyz_path, xyzin_path, gjf_path = write_bdpcs3_outputs(
        working_dir=tmp_path,
        xyzin_path=xyzin,
        result=result,
        bdpcs3_version="legacy",
        basic_charge=0,
        basic_mult=1,
    )

    assert report_path.exists()
    assert xyz_path.exists()
    assert xyzin_path.exists()
    assert gjf_path.exists()

    report_text = report_path.read_text()
    assert "DPCS3 to BDPCS3 report" in report_text
    assert "Rotational constants (MHz)" in report_text

    gjf_text = gjf_path.read_text()
    assert "#HF geom=modredundant Pop=CM5 iop(6/79=1,6/80=1) output=Pickett" in gjf_text
    assert "B 1 2" in gjf_text


def test_bdpcs3_workflow_writes_hbond_target(tmp_path):
    xyzin = tmp_path / "xyzin"
    _write_hbond_xyzin(xyzin)

    result = compute_bdpcs3(xyzin, "updated")
    assert len(result.hbond_targets) == 1
    donor, h_atom, acceptor, r0, r_target, angle_deg, delta, weight = result.hbond_targets[0]
    assert (donor, h_atom, acceptor) == (0, 1, 2)
    assert r0 == pytest.approx(1.84, abs=2.0e-5)
    assert r_target == pytest.approx(1.785, abs=2.0e-5)
    assert angle_deg == pytest.approx(170.0, abs=1.0e-3)
    assert delta < -0.054
    assert weight == pytest.approx(100.0)

    report_path, _xyz_path, _xyzin_path, gjf_path = write_bdpcs3_outputs(
        working_dir=tmp_path,
        xyzin_path=xyzin,
        result=result,
        bdpcs3_version="updated",
        basic_charge=0,
        basic_mult=1,
    )
    assert "Hydrogen-bond targets" in report_path.read_text()
    assert "B 2 3" in gjf_path.read_text()
