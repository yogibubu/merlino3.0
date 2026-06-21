import re
from pathlib import Path

import pytest

from merlino_core.parameters.bdpcs3 import load_bdpcs3_parameters
from merlino_core.parameters.generate_fortran_includes import build_bdpcs3_hbond_include
from survibfit.modify_geom import (
    bdpcs3_hbond_delta,
    bdpcs3_delta_and_order,
    bdpcs3_delta_and_order_updated,
    bdpcs3_function,
    bdpcs3_metric_weights,
    detect_hydrogen_bonds,
)
from survibfit.primitives import Primitive


@pytest.mark.parametrize(
    ("z1", "z2", "r_ang", "expected_delta", "expected_order"),
    [
        (8, 1, 0.95, -0.001056, 1.033895113513574),
        (6, 8, 1.4, -0.0026483056847728136, 0.9672161004820066),
        (9, 6, 1.5, -0.005144190898479566, 0.6065306597126336),
    ],
)
def test_bdpcs3_delta_and_order_matches_reference(
    z1, z2, r_ang, expected_delta, expected_order
):
    delta, bond_order = bdpcs3_delta_and_order(z1, z2, r_ang)
    assert delta == pytest.approx(expected_delta, rel=1e-12)
    assert bond_order == pytest.approx(expected_order, rel=1e-12)


def test_bdpcs3_no_correction_for_weak_cc_cs_bonds():
    delta, bond_order = bdpcs3_delta_and_order(6, 6, 3.0)
    assert bond_order < 0.3
    assert delta == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize(
    ("z1", "z2", "r_ang", "expected_delta"),
    [
        (6, 1, 1.07, -0.0025),  # CH
        (6, 6, 1.52, -0.0024393779762205857),  # CC
        (6, 7, 1.40, -0.0025),  # CN
        (6, 8, 1.42, -0.0025),  # CO
        (6, 16, 1.81, -0.004984244442007778),  # CS
    ],
)
def test_bdpcs3_updated_matches_reference(z1, z2, r_ang, expected_delta):
    delta, _ = bdpcs3_delta_and_order_updated(z1, z2, r_ang)
    assert delta == pytest.approx(expected_delta, abs=5e-6)


def test_bdpcs3_selector_uses_updated_by_default():
    assert bdpcs3_function("updated") is bdpcs3_delta_and_order_updated
    assert bdpcs3_function("unified") is bdpcs3_delta_and_order_updated


def test_bdpcs3_updated_differs_from_legacy():
    legacy, _ = bdpcs3_delta_and_order(6, 16, 1.81)
    updated, _ = bdpcs3_delta_and_order_updated(6, 16, 1.81)
    assert abs(updated - legacy) > 1e-4


def test_bdpcs3_updated_differs_from_legacy_cn():
    legacy, _ = bdpcs3_delta_and_order(6, 7, 1.40)
    updated, _ = bdpcs3_delta_and_order_updated(6, 7, 1.40)
    assert abs(updated - legacy) > 1e-4


def test_bdpcs3_override_does_not_artificially_lower_bond_order():
    # Topology BO can be spuriously low for strained but real sigma bonds.
    # The geometric BO floor prevents disabling corrections in those cases.
    delta_ref, bo_ref = bdpcs3_delta_and_order(6, 6, 1.52)
    delta_override, bo_override = bdpcs3_delta_and_order(
        6, 6, 1.52, bond_order_override=0.273
    )
    assert bo_override == pytest.approx(bo_ref, rel=1e-12)
    assert delta_override == pytest.approx(delta_ref, rel=1e-12)


def test_bdpcs3_ch_not_controlled_by_bond_order_term():
    delta_ref, _ = bdpcs3_delta_and_order(6, 1, 1.11)
    delta_low_bo, _ = bdpcs3_delta_and_order(6, 1, 1.11, bond_order_override=0.05)
    assert delta_low_bo == pytest.approx(delta_ref, rel=1e-12)


def test_bdpcs3_hbond_delta_has_angle_gate_and_distance_damping():
    short_delta = bdpcs3_hbond_delta(8, 8, 1.85, 175.0)
    assert short_delta == pytest.approx(-0.055, abs=1e-8)
    assert bdpcs3_hbond_delta(8, 8, 1.85, 130.0) == pytest.approx(0.0)
    assert abs(bdpcs3_hbond_delta(8, 8, 3.8, 175.0)) < 1.0e-8
    assert bdpcs3_hbond_delta(7, 8, 1.85, 175.0) == pytest.approx(short_delta)


def test_detect_hydrogen_bond_requires_directional_xhy_angle():
    Z = [8, 1, 8]
    covalent = [(0, 1)]
    linear = [
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [2.80, 0.0, 0.0],
    ]
    hbonds = detect_hydrogen_bonds(Z, linear, covalent)
    assert [(hb.donor, hb.hydrogen, hb.acceptor) for hb in hbonds] == [(0, 1, 2)]
    assert hbonds[0].angle_deg == pytest.approx(180.0)

    bent = [
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [0.96, 1.84, 0.0],
    ]
    assert detect_hydrogen_bonds(Z, bent, covalent) == []


def test_bdpcs3_metric_weights_distinguish_hbond_from_stretch():
    prims = [
        Primitive("bond", (0, 1)),
        Primitive("bond", (1, 2)),
        Primitive("angle", (0, 1, 2)),
    ]
    weights = bdpcs3_metric_weights(
        prims,
        [8, 1, 8],
        [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [2.80, 0.0, 0.0]],
        hbond_pairs={(1, 2)},
    )
    assert weights.tolist() == pytest.approx([1000.0, 100.0, 100.0])


def test_bdpcs3_fortran_include_matches_shared_parameters():
    params = load_bdpcs3_parameters()
    root = Path(__file__).resolve().parents[2]
    include_path = root / "merlino_core" / "parameters" / "fortran" / "bdpcs3_hbond_params.inc"
    text = include_path.read_text()
    assert text == build_bdpcs3_hbond_include()
    values = {}
    for name, value in re.findall(r"PARAMETER\s*\(([^=]+)=([^)]+)\)", text):
        values[name.strip()] = float(value.strip().replace("D", "E"))

    assert values["BDPCS3_HB_ANGLE_MIN"] == pytest.approx(params.hbond.angle_threshold_deg)
    assert values["BDPCS3_HB_DIST_CUTOFF"] == pytest.approx(params.hbond.distance_cutoff_ang)
    assert values["BDPCS3_HB_DIST_WIDTH"] == pytest.approx(params.hbond.distance_width_ang)
    assert values["BDPCS3_HB_SEARCH_CUTOFF"] == pytest.approx(params.hbond.search_cutoff_ang)
    assert values["BDPCS3_HB_DELTA_OO"] == pytest.approx(params.hbond.correction_ang(8, 8))
    assert values["BDPCS3_HB_DELTA_NN"] == pytest.approx(params.hbond.correction_ang(7, 7))
    assert values["BDPCS3_W_STRETCH"] == pytest.approx(params.weights.stretch)
    assert values["BDPCS3_W_ANGLE"] == pytest.approx(params.weights.angle)
    assert values["BDPCS3_W_HBOND"] == pytest.approx(params.weights.hbond)
    assert values["BDPCS3_W_TORSION_MIN"] == pytest.approx(params.weights.torsion_min)
    assert values["BDPCS3_W_FRAGMENT"] == pytest.approx(params.weights.fragment)
