from __future__ import annotations

import pytest
import numpy as np

from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive
from merlino_semiexp import (
    CorrectedRotationalConstants,
    IsotopologueObservation,
    QMParameterPredicate,
    RotationalConstants,
    SemiexperimentalFitRequest,
    VibrationalCorrection,
    corrected_constants_rows,
    fit_semiexperimental_geometry,
    parse_substitutions,
    read_observations_csv,
    write_observations_csv,
)
from merlino_core import repo_root
from merlino_vpt2_vci import (
    DavidsonSettings,
    ForceFieldSource,
    VCIRequest,
    inventory_vpt2_vci_backends,
)


def test_davidson_settings_validation():
    settings = DavidsonSettings(n_roots=3, max_subspace=12, max_iter=50, convergence=1.0e-7)
    settings.validate()

    with pytest.raises(ValueError):
        DavidsonSettings(n_roots=5, max_subspace=4).validate()


def test_vci_request_validation(tmp_path):
    request = VCIRequest(
        force_field=ForceFieldSource(tmp_path / "qff.log"),
        max_quanta=4,
        basis_energy_cutoff_cm=8000.0,
    )
    request.validate()

    with pytest.raises(ValueError):
        VCIRequest(force_field=ForceFieldSource(tmp_path / "qff.log"), max_quanta=0).validate()


def test_vpt2_vci_inventory_records_active_backend_status():
    inventory = inventory_vpt2_vci_backends(repo_root(__file__))
    assert inventory.harmonic_source is not None
    assert inventory.harmonic_source.name == "gf_core.f"
    assert {path.name for path in inventory.active_fortran_sources} == {
        "davidson_core.f",
        "gf_core.f",
        "vci_core.f",
        "vpt2_core.f",
    }
    assert inventory.davidson_backend is not None
    assert inventory.davidson_backend.name == "davidson_core.f"
    notes = " ".join(inventory.notes)
    assert "GF" in notes
    assert "Davidson" in notes

def test_semiexperimental_correction_subtracts_vibrational_delta():
    observed = RotationalConstants(1000.0, 800.0, 600.0)
    correction = VibrationalCorrection(1.0, -2.0, 0.5, source="qm")
    corrected = CorrectedRotationalConstants(observed, correction).equilibrium

    assert corrected.as_tuple() == (999.0, 802.0, 599.5)


def test_semiexperimental_fit_request_validation(tmp_path):
    obs = IsotopologueObservation(
        label="parent",
        constants=RotationalConstants(1000.0, 800.0, 600.0),
    )
    request = SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs,))
    request.validate()

    duplicate = SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs, obs))
    with pytest.raises(ValueError):
        duplicate.validate()


def test_semiexperimental_observations_csv_roundtrip(tmp_path):
    observations = (
        IsotopologueObservation(
            label="parent",
            constants=RotationalConstants(1000.0, 800.0, 600.0),
            correction=VibrationalCorrection(1.0, 2.0, 3.0, source="gaussian"),
            weights=RotationalConstants(100.0, 25.0, 4.0),
        ),
        IsotopologueObservation(
            label="13C1",
            constants=RotationalConstants(990.0, 790.0, 590.0),
            substitutions={1: 13},
            correction=VibrationalCorrection(0.5, 1.5, 2.5, source="gaussian"),
        ),
    )

    path = write_observations_csv(tmp_path / "observations.csv", observations)
    loaded = read_observations_csv(path)
    rows = corrected_constants_rows(loaded)

    assert loaded[1].substitutions == {1: 13}
    assert loaded[0].weights is not None
    assert loaded[0].weights.as_tuple() == pytest.approx((100.0, 25.0, 4.0))
    assert rows[0]["A_e_MHz"] == 999.0
    assert rows[1]["C_e_MHz"] == 587.5


def test_semiexperimental_substitution_parser():
    assert parse_substitutions("2:13;5:18") == {2: 13, 5: 18}
    with pytest.raises(ValueError):
        parse_substitutions("0:13")


def test_analytic_b_matrix_for_bond_primitive():
    coords = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=float)
    b = b_matrix_analytic([Primitive("bond", (0, 1))], coords)

    assert b.shape == (1, 6)
    assert b[0] == pytest.approx([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0])


def test_semiexperimental_geometry_fit_reduces_rotational_residuals(tmp_path):
    atoms = ["O", "H", "H"]
    target = np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 0.9572],
            [0.9266, 0.0000, -0.2396],
        ],
        dtype=float,
    )
    initial = target.copy()
    initial[1, 2] += 0.08
    initial[2, 0] -= 0.05
    xyz = tmp_path / "water_initial.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "distorted water",
                *[f"{atom} {x:.8f} {y:.8f} {z:.8f}" for atom, (x, y, z) in zip(atoms, initial)],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    parent_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, target)))
    d1_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, target, [None, 2, None])))
    observations = (
        IsotopologueObservation("parent", parent_constants),
        IsotopologueObservation("D1", d1_constants, substitutions={2: 2}),
    )
    request = SemiexperimentalFitRequest(xyz, observations)
    initial_rms = _rotconst_rms(atoms, initial, observations)

    result = fit_semiexperimental_geometry(request, max_iter=8, outdir=tmp_path / "semiexp")

    assert result.rms_MHz < initial_rms
    assert result.diagnostics.observable == "moments"
    assert result.diagnostics.components == ("Ia", "Ib", "Ic")
    assert result.b_matrix.shape[0] == len(result.gic_labels)
    assert result.b_matrix.shape[1] == 3 * len(atoms)
    assert result.hessian.shape == result.covariance.shape
    assert result.correlation.shape == result.covariance.shape
    assert result.stationary_point in {"minimum", "flat_or_rank_deficient"}
    assert result.diagnostics.rank <= result.jacobian.shape[1]
    assert result.diagnostics.accepted_steps >= 0
    assert result.diagnostics.rejected_steps >= 0
    assert result.diagnostics.convergence_reason in {
        "rms_tolerance",
        "gradient_tolerance",
        "objective_tolerance",
        "max_iter",
    }
    assert all(np.isfinite(parameter.sigma) for parameter in result.parameters)
    assert (tmp_path / "semiexp" / "semiexp_geometry.xyz").exists()
    assert (tmp_path / "semiexp" / "semiexp_parameters.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_residuals.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_covariance.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_correlation.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_hessian.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_hessian_eigenvalues.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_diagnostics.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_manifest.json").exists()


def test_semiexperimental_qm_predicate_adds_weighted_parameter_prior(tmp_path):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9572], [0.9266, 0.0, -0.2396]])
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords))),
    )
    request = SemiexperimentalFitRequest(
        xyz,
        (observation,),
        qm_predicates=(QMParameterPredicate("GIC001", 1.0, 0.1, source="qm-estimate"),),
    )

    result = fit_semiexperimental_geometry(request, max_iter=1)

    assert any(residual.isotopologue == "qm-estimate" for residual in result.residuals)
    assert result.diagnostics.rank <= result.jacobian.shape[1]


def test_planar_rotational_constants_auto_selects_stable_pair(tmp_path):
    atoms = ["C", "O", "H", "H"]
    coords = np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [1.2000, 0.0000, 0.0000],
            [-0.6000, 0.9000, 0.0000],
            [-0.6000, -0.9000, 0.0000],
        ],
        dtype=float,
    )
    xyz = tmp_path / "formaldehyde.xyz"
    xyz.write_text(
        "\n".join(["4", "planar", *[f"{a} {x:.8f} {y:.8f} {z:.8f}" for a, (x, y, z) in zip(atoms, coords)]])
        + "\n",
        encoding="utf-8",
    )
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords))),
    )

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(
            xyz,
            (observation,),
            observable="rotational_constants",
            rotational_components="auto",
        ),
        max_iter=1,
    )

    assert result.diagnostics.planar is True
    assert len(result.diagnostics.components) == 2
    assert set(result.diagnostics.components).issubset({"A", "B", "C"})


def test_semiexperimental_fit_honors_fixed_gic_parameters(tmp_path):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9572], [0.9266, 0.0, -0.2396]])
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords))),
    )

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(xyz, (observation,), fixed_parameters=("GIC001",)),
        max_iter=1,
    )

    assert result.parameters[0].active is False
    assert any(parameter.active for parameter in result.parameters)


def _structure(atoms, coords, isotopes=None):
    return Structure.from_atoms_coords(atoms, [tuple(row) for row in coords], isotopes=isotopes)


def _rotconst_rms(atoms, coords, observations):
    residuals = []
    for observation in observations:
        isotopes = [None] * len(atoms)
        for atom_index, mass_number in observation.substitutions.items():
            isotopes[atom_index - 1] = mass_number
        calculated = rotational_constants_MHz(_structure(atoms, coords, isotopes))
        residuals.extend(obs - calc for obs, calc in zip(observation.corrected.as_tuple(), calculated))
    return float(np.sqrt(np.mean(np.square(residuals))))
