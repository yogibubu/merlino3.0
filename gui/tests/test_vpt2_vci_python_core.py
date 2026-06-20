from __future__ import annotations

import numpy as np

from merlino_vpt2_vci import (
    QuarticForceField,
    davidson_lowest,
    generate_vibrational_basis,
    lower_to_symmetric,
    read_gaussian_fchk_qff,
    read_indexed_qff_text,
    run_python_vci_from_gaussian_fchk,
    solve_vci,
    solve_wilson_gf,
    zero_anharmonic_force_field,
)


def test_wilson_gf_identity_case_returns_square_root_frequencies():
    result = solve_wilson_gf(np.diag([4.0, 9.0]), np.eye(2))

    assert np.allclose(result.eigenvalues, [4.0, 9.0])
    assert np.allclose(result.frequencies_cm, [2.0, 3.0])


def test_vci_harmonic_basis_and_energies_are_deterministic():
    qff = zero_anharmonic_force_field(np.array([100.0, 200.0]))
    result = solve_vci(qff, max_quanta=1)

    assert result.basis == ((0, 0), (0, 1), (1, 0))
    assert np.allclose(result.energies_cm, [150.0, 250.0, 350.0])


def test_vci_quartic_term_shifts_oscillator_ground_state():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([100.0]),
        cubic_cm={},
        quartic_cm={(0, 0, 0, 0): 4.0},
    )
    result = solve_vci(qff, max_quanta=0)

    assert np.allclose(result.energies_cm, [53.0])


def test_gaussian_fchk_qff_reader_uses_real_anharmonic_blocks():
    data = read_gaussian_fchk_qff(__import__("pathlib").Path("gui/tests/gaussian/h2o.fchk"))

    assert data.masses_amu.shape == (3,)
    assert data.cartesian_hessian_lower.shape == (45,)
    assert np.allclose(data.anharmonic_frequencies_cm[:3], [2123.50470, 4016.61987, 4266.73074])
    assert lower_to_symmetric(data.cartesian_hessian_lower).shape == (9, 9)


def test_python_fchk_workflow_produces_vci_levels():
    run = run_python_vci_from_gaussian_fchk(__import__("pathlib").Path("gui/tests/gaussian/h2o.fchk"), max_quanta=1)

    assert run.gf.frequencies_cm.size == 9
    assert len(run.vci.basis) == len(generate_vibrational_basis(3, 1))
    assert run.vci.energies_cm[0] > 0.0


def test_indexed_qff_text_feeds_quartic_vci(tmp_path):
    qff_file = tmp_path / "field.qff"
    qff_file.write_text(
        "\n".join(
            [
                "FREQ 1 100.0",
                "QUARTIC 1 1 1 1 4.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    qff = read_indexed_qff_text(qff_file)
    result = solve_vci(qff, max_quanta=0)

    assert np.allclose(result.energies_cm, [53.0])


def test_independent_davidson_matches_dense_symmetric_diagonalization():
    matrix = np.array(
        [
            [2.0, 0.1, 0.0, 0.0],
            [0.1, 3.0, 0.2, 0.0],
            [0.0, 0.2, 5.0, 0.3],
            [0.0, 0.0, 0.3, 8.0],
        ]
    )
    expected, _ = np.linalg.eigh(matrix)

    result = davidson_lowest(lambda vector: matrix @ vector, np.diag(matrix), n_roots=2)

    assert result.converged
    assert np.allclose(result.eigenvalues, expected[:2], atol=1.0e-8)


def test_vci_davidson_path_matches_dense_path():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([100.0, 140.0]),
        cubic_cm={(0, 0, 1): 1.5},
        quartic_cm={(0, 0, 0, 0): 0.2, (0, 1, 1, 1): -0.1},
    )

    dense = solve_vci(qff, max_quanta=3, n_roots=3)
    iterative = solve_vci(qff, max_quanta=3, n_roots=3, method="davidson")

    assert iterative.davidson is not None
    assert iterative.davidson.converged
    assert np.allclose(iterative.energies_cm, dense.energies_cm, atol=1.0e-7)
