from __future__ import annotations

import numpy as np
import pytest

from merlino_fit.survibfit.primitives import Primitive
from merlino_gf import (
    BOHR_TO_ANGSTROM,
    gf_from_hessian_input_and_gic_definition,
    gf_from_gaussian_fchk_with_merlino_gics,
    gf_from_hessian_input_with_merlino_gics,
    pulay_scale_internal_hessian,
    run_gic_gf_report_from_fchk,
    run_gf_report_from_fchk,
    solve_wilson_gf,
)
from merlino_gic import GICDefinition
from merlino_vpt2_vci.gaussian_qff import hessian_input_from_gaussian_fchk
from merlino_vpt2_vci import (
    QuarticForceField,
    AnharmonicInput,
    VCIOptions,
    ScientificValidationError,
    compare_vpt2_vci,
    davidson_lowest,
    generate_vibrational_basis,
    lower_to_symmetric,
    load_force_field,
    read_gaussian_fchk_qff,
    read_indexed_qff_text,
    run_python_vci_from_gaussian_fchk,
    run_vpt2_vci_report,
    solve_vci,
    solve_vci_from_anharmonic_input,
    solve_vpt2_from_anharmonic_input,
    validate_force_field,
    zero_anharmonic_force_field,
)


def test_wilson_gf_identity_case_returns_square_root_frequencies():
    result = solve_wilson_gf(np.diag([4.0, 9.0]), np.eye(2))

    assert np.allclose(result.eigenvalues, [4.0, 9.0])
    assert np.allclose(result.frequencies_cm, [2.0, 3.0])


def test_vci_harmonic_basis_and_energies_are_deterministic():
    qff = zero_anharmonic_force_field(np.array([100.0, 200.0]))
    result = solve_vci(qff, max_quanta=1)

    assert result.basis == ((0, 0), (1, 0), (0, 1))
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

    assert data.atomic_numbers.tolist() == [1, 8, 1]
    assert data.cartesian_coordinates_bohr.shape == (3, 3)
    assert data.masses_amu.shape == (3,)
    assert data.cartesian_hessian_lower.shape == (45,)
    assert np.allclose(data.anharmonic_frequencies_cm[:3], [2123.50470, 4016.61987, 4266.73074])
    assert lower_to_symmetric(data.cartesian_hessian_lower).shape == (9, 9)
    hessian_input = data.to_hessian_input()
    assert hessian_input.source == "gaussian-fchk"
    assert hessian_input.cartesian_hessian.shape == (9, 9)


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


def test_canonical_anharmonic_input_runs_vci_and_shifts_levels():
    harmonic = solve_vci_from_anharmonic_input(
        AnharmonicInput(
            harmonic_frequencies_cm=np.array([1000.0, 1500.0]),
            anharmonic_frequencies_cm=np.array([]),
            cubic_cm={},
            quartic_cm={},
        ),
        max_quanta=2,
        n_roots=4,
    )
    anharmonic = solve_vci_from_anharmonic_input(
        AnharmonicInput(
            harmonic_frequencies_cm=np.array([1000.0, 1500.0]),
            anharmonic_frequencies_cm=np.array([]),
            cubic_cm={(0, 0, 1): -25.0},
            quartic_cm={(0, 0, 0, 0): 8.0, (0, 0, 1, 1): -3.0},
        ),
        max_quanta=2,
        n_roots=4,
    )

    assert np.allclose(harmonic.excitation_energies_cm[:4], [0.0, 1000.0, 1500.0, 2000.0])
    assert not np.allclose(anharmonic.excitation_energies_cm, harmonic.excitation_energies_cm)
    assert np.all(anharmonic.excitation_energies_cm[1:] > 0.0)


def test_vci_basis_cutoff_frequency_window_pruning_blocks_and_contributions():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([500.0, 1000.0, 1800.0]),
        cubic_cm={(0, 0, 1): 0.001, (1, 1, 1): 20.0},
        quartic_cm={(1, 1, 1, 1): 5.0, (0, 0, 1, 1): 0.0001},
    )
    options = VCIOptions(
        frequency_min_cm=800.0,
        frequency_max_cm=1600.0,
        basis_energy_cutoff_cm=2200.0,
        force_constant_threshold_cm=0.01,
        mode_symmetries=("a", "b", "c"),
        separate_symmetry_blocks=True,
        coefficient_threshold=0.05,
    )

    result = solve_vci(qff, max_quanta=4, n_roots=4, options=options)

    assert len(result.basis) == 3
    assert result.basis == ((0,), (1,), (2,))
    assert {block.label for block in result.blocks} == {"A", "b"}
    assert len(result.state_contributions) == 3
    assert np.allclose(result.state_contributions[1].mode_quanta, [1.0], atol=1.0e-8)
    assert result.state_contributions[1].dominant_basis_states[0][0] == (1,)


def test_vci_per_mode_and_excitation_class_limits():
    basis = generate_vibrational_basis(
        3,
        max_quanta=4,
        frequencies_cm=np.array([100.0, 200.0, 300.0]),
        mode_max_quanta=(3, 2, 1),
        excitation_class_limits={
            1: (1, 2),
            2: (2, 3),
            3: (3, 3),
        },
    )

    assert (0, 0, 0) in basis
    assert (3, 0, 0) not in basis
    assert (2, 2, 0) not in basis
    assert (2, 1, 0) in basis
    assert (1, 1, 1) in basis
    assert (2, 1, 1) not in basis
    assert all(state[1] <= 2 and state[2] <= 1 for state in basis)


def test_vpt2_quartic_first_order_matches_ground_state_shift():
    result = solve_vpt2_from_anharmonic_input(
        AnharmonicInput(
            harmonic_frequencies_cm=np.array([100.0]),
            anharmonic_frequencies_cm=np.array([]),
            cubic_cm={},
            quartic_cm={(0, 0, 0, 0): 4.0},
        ),
        max_quanta=0,
    )

    assert np.allclose(result.energies_cm, [53.0])
    assert np.allclose(result.states[0].first_order_cm, 3.0)


def test_vpt2_vci_comparison_uses_same_reduced_mode_selection():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([500.0, 1000.0, 1500.0]),
        cubic_cm={(1, 1, 2): -2.0, (0, 1, 1): 100.0},
        quartic_cm={(1, 1, 1, 1): 0.8, (2, 2, 2, 2): 0.5, (0, 0, 1, 1): 200.0},
    )
    options = VCIOptions(
        active_modes=(1, 2),
        mode_max_quanta=(0, 3, 2),
        excitation_class_limits={1: (1, 2), 2: (2, 3)},
        force_constant_threshold_cm=1.0,
    )

    comparison = compare_vpt2_vci(qff, max_quanta=3, n_roots=4, options=options)

    assert comparison.vpt2.basis == comparison.vci.basis[: len(comparison.vpt2.basis)]
    assert all(len(state) == 2 for state in comparison.vpt2.basis)
    assert all(state[0] <= 3 and state[1] <= 2 for state in comparison.vpt2.basis)
    assert (0, 0) in comparison.vpt2.basis
    assert comparison.energy_differences_cm.shape == (4,)
    assert np.max(np.abs(comparison.excitation_differences_cm)) < 20.0


def test_gf_from_gaussian_cartesian_hessian_uses_merlino_nonredundant_gics():
    path = __import__("pathlib").Path("gui/tests/gaussian/h2o.fchk")
    canonical = hessian_input_from_gaussian_fchk(path)
    result = gf_from_hessian_input_with_merlino_gics(canonical)
    adapter_result = gf_from_gaussian_fchk_with_merlino_gics(path)

    assert result.b_matrix.shape == (3, 9)
    assert result.force_constants.shape == (3, 3)
    assert result.g_matrix.shape == (3, 3)
    assert result.ped.values.shape == (3, 3)
    assert "bond(1,2)" in result.gic_labels[0]
    assert "angle(1,2,3)" in result.gic_labels[2]
    assert np.all(result.frequencies_cm > 0.0)
    assert np.allclose(result.ped.values.sum(axis=0), np.full(3, 100.0))
    assert np.allclose(result.frequencies_cm, [2169.878, 4141.256, 4392.363], atol=1.0e-3)
    assert np.allclose(adapter_result.frequencies_cm, result.frequencies_cm)


def test_frozen_gic_gf_branch_scales_internal_hessian_pulay_style(tmp_path):
    path = __import__("pathlib").Path("gui/tests/gaussian/h2o.fchk")
    canonical = hessian_input_from_gaussian_fchk(path)
    definition = GICDefinition(
        atom_symbols=("H", "O", "H"),
        atomic_numbers=(1, 8, 1),
        reference_coordinates_angstrom=tuple(tuple(row) for row in canonical.cartesian_coordinates_bohr * BOHR_TO_ANGSTROM),
        primitives=(
            Primitive("bond", (0, 1)),
            Primitive("bond", (1, 2)),
            Primitive("angle", (0, 1, 2)),
        ),
        u_matrix=np.eye(3),
        labels=("GIC001 R(1,2)", "GIC002 R(2,3)", "GIC003 A(1,2,3)"),
        names=("R1", "R2", "A1"),
        irreps=("A1", "A1", "A1"),
        point_group="C2v",
    )

    unscaled = gf_from_hessian_input_and_gic_definition(canonical, definition)
    factors = np.array([1.0, 0.95, 0.9])
    scaled = gf_from_hessian_input_and_gic_definition(canonical, definition, scaling_factors=factors)
    expected = unscaled.force_constants * np.sqrt(np.outer(factors, factors))

    assert unscaled.b_matrix.shape == (3, 9)
    assert np.allclose(scaled.force_constants, expected)
    assert np.allclose(pulay_scale_internal_hessian(unscaled.force_constants, factors), expected)
    assert not np.allclose(scaled.frequencies_cm, unscaled.frequencies_cm)
    assert np.allclose(scaled.ped.values.sum(axis=0), np.full(3, 100.0))

    schema = definition.write(tmp_path / "gic_definition.json")
    scale = tmp_path / "scale.txt"
    scale.write_text("GIC003 0.90\nR2 0.95\n", encoding="utf-8")
    report = run_gic_gf_report_from_fchk(path, schema, scale_path=scale)

    assert "Frozen GIC definition" in report.text
    assert "Pulay Hessian scaling: applied" in report.text
    assert report.result.scaling_factors.tolist() == pytest.approx([1.0, 0.95, 0.9])


def test_gui_service_reports_are_independent_from_qt(tmp_path):
    fchk = __import__("pathlib").Path("gui/tests/gaussian/h2o.fchk")
    gf_report = run_gf_report_from_fchk(fchk)

    assert "GF/PED from Merlino non-redundant GICs" in gf_report.text
    assert gf_report.result.ped.values.shape == (3, 3)

    qff_file = tmp_path / "field.qff"
    qff_file.write_text(
        "\n".join(
            [
                "FREQ 1 1000.0",
                "FREQ 2 1500.0",
                "QUARTIC 1 1 1 1 0.8",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    qff = load_force_field(qff_path=qff_file)
    report = run_vpt2_vci_report(qff, max_quanta=2, roots=3, options=VCIOptions(active_modes=(0, 1)))

    assert "VPT2/VCI comparison" in report.text
    assert len(report.comparison.vci.basis) >= 3


def test_force_field_validation_rejects_bad_modes():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([1000.0]),
        cubic_cm={(0, 0, 2): 1.0},
        quartic_cm={},
    )

    with pytest.raises(ScientificValidationError):
        validate_force_field(qff)
