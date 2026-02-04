import pytest

from survibfit.modify_geom import bdpcs3_delta_and_order, bdpcs3_delta_and_order_updated


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
        (6, 1, 1.07, -0.002078154),  # CH
        (6, 6, 1.52, -0.002951739360067995),  # CC
        (6, 8, 1.42, -0.0017679),  # CO
        (6, 16, 1.81, -0.006008788170109689),  # CS
    ],
)
def test_bdpcs3_updated_matches_reference(z1, z2, r_ang, expected_delta):
    delta, _ = bdpcs3_delta_and_order_updated(z1, z2, r_ang)
    assert delta == pytest.approx(expected_delta, abs=5e-6)


def test_bdpcs3_updated_differs_from_legacy():
    legacy, _ = bdpcs3_delta_and_order(6, 16, 1.81)
    updated, _ = bdpcs3_delta_and_order_updated(6, 16, 1.81)
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
