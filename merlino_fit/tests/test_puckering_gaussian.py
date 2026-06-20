import numpy as np

from survibfit.puckering_gaussian import (
    angular_step_to_target,
    auto_ring_indices,
    build_gjf_links,
    canonical_ring_indices,
    four_ring_target_gic,
    parse_ring_indices,
    puckering_state,
    ring_puckering_gic_lines,
)


def test_parse_ring_indices_validates():
    assert parse_ring_indices("1,2,3,4,5", 8) == [0, 1, 2, 3, 4]


def test_ring_numbering_is_canonical_under_rotation_and_reversal():
    assert canonical_ring_indices([2, 3, 4, 0, 1]) == [0, 1, 2, 3, 4]
    assert canonical_ring_indices([3, 2, 1, 0, 4]) == [0, 1, 2, 3, 4]
    assert ring_puckering_gic_lines([2, 3, 4, 0, 1]) == ring_puckering_gic_lines([3, 2, 1, 0, 4])


def test_five_ring_gic_contains_q_and_phi():
    coords = np.array(
        [
            [1.0, 0.0, 0.15],
            [0.30901699, 0.95105652, -0.10],
            [-0.80901699, 0.58778525, 0.12],
            [-0.80901699, -0.58778525, -0.08],
            [0.30901699, -0.95105652, -0.09],
        ],
        dtype=float,
    )
    atoms = ["C"] * 5
    lines, manifest = build_gjf_links(
        atoms,
        coords,
        [0, 1, 2, 3, 4],
        0.0,
        20.0,
        10.0,
        chk_prefix="test_phi",
    )
    text = "\n".join(lines)
    assert "RPck001(Inactive)=" in text
    assert "RPck002(Inactive)=" in text
    assert "QPck001=SQRT(RPck001*RPck001+RPck002*RPck002)" in text
    assert "PhiP001(NSteps=1,StepSize=" in text
    assert "--Link1--" in text
    assert len(manifest) == 3


def test_auto_ring_indices_detects_five_membered_ring():
    coords = np.array(
        [
            [1.20630072, 0.18231452, 0.18],
            [0.19937601, 1.20359844, -0.12],
            [-1.08307957, 0.56155022, 0.15],
            [-0.86875599, -0.85654132, -0.10],
            [0.54615884, -1.09092187, -0.11],
        ],
        dtype=float,
    )
    atoms = ["O", "C", "C", "C", "C"]
    assert auto_ring_indices(atoms, coords) == [0, 1, 2, 3, 4]


def test_four_ring_gic_contains_q_and_phi():
    lines, step_deg = four_ring_target_gic([0, 1, 2, 3], 15.0, -5.0)
    text = "\n".join(lines)
    assert "RPck001(Inactive)=D(2,3,4,1)" in text
    assert "QPck001=SQRT(RPck001*RPck001)" in text
    assert "PhiP001(NSteps=1,StepSize=" in text
    assert step_deg == 20.0


def test_seven_ring_emits_two_puckering_pairs():
    theta = np.linspace(0.0, 2.0 * np.pi, 7, endpoint=False)
    coords = np.column_stack(
        [
            np.cos(theta),
            np.sin(theta),
            0.1 * np.sin(2.0 * theta),
        ]
    )
    atoms = ["C"] * 7
    lines, manifest = build_gjf_links(
        atoms,
        coords,
        list(range(7)),
        0.0,
        0.0,
        10.0,
        chk_prefix="test_phi7",
    )
    text = "\n".join(lines)
    assert "T007(Inactive)=D(7,1,2,3)" in text
    assert "RPck004(Inactive)=" in text
    assert "QPck001=SQRT(RPck001*RPck001+RPck002*RPck002)" in text
    assert "PhiP001(NSteps=1,StepSize=" in text
    assert "QPck002=SQRT(RPck003*RPck003+RPck004*RPck004)" in text
    assert "PhiP002=ATAN2(RPck004,RPck003)" in text
    assert len(manifest) == 1


def test_puckering_state_and_angular_step():
    coords = np.array(
        [
            [1.0, 0.0, 0.15],
            [0.30901699, 0.95105652, -0.10],
            [-0.80901699, 0.58778525, 0.12],
            [-0.80901699, -0.58778525, -0.08],
            [0.30901699, -0.95105652, -0.09],
        ],
        dtype=float,
    )
    state = puckering_state(coords, [0, 1, 2, 3, 4])
    assert state.q > 0.0
    assert 0.0 <= state.phi_deg < 360.0
    assert angular_step_to_target(10.0, 350.0) == 20.0
