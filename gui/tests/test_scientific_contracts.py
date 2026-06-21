from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np

from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive
from merlino_gf import BOHR_TO_ANGSTROM
from merlino_semiexp import (
    CorrectedRotationalConstants,
    DEFAULT_SEMIEXP_OBSERVABLE,
    DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
    ElectronicCorrection,
    HYDROGEN_PARAMETER_CONSTRAINT,
    IsotopologueObservation,
    ParameterClassConstraint,
    QMParameterPredicate,
    RotationalConstants,
    SEMIEXP_JOB_SCHEMA,
    SemiexperimentalFitRequest,
    VibrationalCorrection,
    cartesian_symmetry_coordinate_model,
    corrected_constants_rows,
    fit_semiexperimental_geometry,
    kraitchman_comparison,
    kraitchman_seed_geometry,
    parse_substitutions,
    preview_semiexperimental_conditioning,
    preview_semiexperimental_gics,
    read_geometry_input,
    read_observations,
    read_observations_csv,
    read_semiexperimental_job,
    semiexperimental_latex_tables,
    validate_semiexperimental_request,
    write_semiexperimental_html_report,
    write_observations_csv,
)
from merlino_vpt2_vci.gaussian_qff import hessian_input_from_gaussian_fchk
from merlino_core import repo_root
from merlino_semiexp.fit import (
    _atomic_number,
    _fixed_primitives_from_patterns,
    _gic_model,
    _hydrogen_fixed_primitives,
    _make_gicforge_backend,
    _primitive_constraint_key,
    _symmetry_expanded_fixed_primitives,
)
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
    electronic = ElectronicCorrection(0.25, 0.5, -0.25, source="electronic")
    corrected = CorrectedRotationalConstants(observed, correction, electronic).equilibrium

    assert corrected.as_tuple() == (998.75, 801.5, 599.75)


def test_semiexperimental_correction_can_use_msr_additive_convention():
    observed = RotationalConstants(1000.0, 800.0, 600.0)
    correction = VibrationalCorrection(1.0, 2.0, 3.0, source="MSR DBvib", convention="additive")
    electronic = ElectronicCorrection(0.25, 0.5, 0.75, source="MSR DBEle", convention="additive")
    corrected = CorrectedRotationalConstants(observed, correction, electronic).equilibrium

    assert corrected.as_tuple() == pytest.approx((1001.25, 802.5, 603.75))


def test_semiexperimental_fit_request_validation(tmp_path):
    obs = IsotopologueObservation(
        label="parent",
        constants=RotationalConstants(1000.0, 800.0, 600.0),
    )
    request = SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs,))
    request.validate()
    assert request.observable == DEFAULT_SEMIEXP_OBSERVABLE == "moments"
    assert request.rotational_components == DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS == "auto"

    duplicate = SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs, obs))
    with pytest.raises(ValueError):
        duplicate.validate()


def test_semiexperimental_observations_csv_roundtrip(tmp_path):
    observations = (
        IsotopologueObservation(
            label="parent",
            constants=RotationalConstants(1000.0, 800.0, 600.0),
            correction=VibrationalCorrection(1.0, 2.0, 3.0, source="gaussian"),
            electronic_correction=ElectronicCorrection(0.1, 0.2, 0.3, source="rel"),
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
    assert rows[0]["A_e_MHz"] == pytest.approx(998.9)
    assert rows[1]["C_e_MHz"] == 587.5


def test_semiexperimental_structured_observations_toml_and_json(tmp_path):
    toml_path = tmp_path / "isotopologues.toml"
    toml_path.write_text(
        """
[[isotopologues]]
label = "parent"
substitutions = ""
[isotopologues.constants]
A_MHz = 1000.0
B_MHz = 800.0
C_MHz = 600.0
[isotopologues.vibrational_correction]
delta_A_MHz = 1.0
delta_B_MHz = 2.0
delta_C_MHz = 3.0
source = "vib"
[isotopologues.electronic_correction]
delta_A_MHz = 0.1
delta_B_MHz = 0.2
delta_C_MHz = 0.3
source = "elec"
[isotopologues.sigma_MHz]
A_MHz = 0.01
B_MHz = 0.02
C_MHz = 0.05

[[isotopologues]]
label = "D2"
substitutions = [{ atom = 2, mass = "D" }]
[isotopologues.constants]
A_MHz = 900.0
B_MHz = 700.0
C_MHz = 500.0
""",
        encoding="utf-8",
    )
    json_path = tmp_path / "isotopologues.json"
    json_path.write_text(
        """
{
  "isotopologues": [
    {
      "label": "13C1",
      "substitutions": {"1": 13},
      "constants": {"A_MHz": 990.0, "B_MHz": 790.0, "C_MHz": 590.0},
      "vibrational_correction": {"delta_A_MHz": 0.5, "delta_B_MHz": 1.5, "delta_C_MHz": 2.5}
    }
  ]
}
""",
        encoding="utf-8",
    )

    toml_obs = read_observations(toml_path)
    json_obs = read_observations(json_path)

    assert toml_obs[0].corrected.as_tuple() == pytest.approx((998.9, 797.8, 596.7))
    assert toml_obs[0].weights is not None
    assert toml_obs[1].substitutions == {2: 2}
    assert json_obs[0].substitutions == {1: 13}
    assert json_obs[0].corrected.C_MHz == pytest.approx(587.5)


def test_semiexperimental_reads_gaussian_cartesian_com_and_modredundant_constraints(tmp_path):
    gaussian_input = tmp_path / "water.com"
    gaussian_input.write_text(
        """
#p hf/sto-3g opt=modredundant

water constrained

0 1
O  0.00000000  0.00000000  0.00000000
H  0.00000000  0.00000000  0.95720000
H  0.92660000  0.00000000 -0.23960000

B 1 2 F
A 2 1 3 F
""",
        encoding="utf-8",
    )

    geometry = read_geometry_input(gaussian_input)

    assert geometry.source_format == "gaussian_cartesian_modredundant"
    assert geometry.comment == "water constrained"
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.coordinates_angstrom.shape == (3, 3)
    assert "bond(1,2)" in geometry.fixed_parameters
    assert "angle(2,1,3)" in geometry.fixed_parameters


def test_semiexperimental_rejects_gaussian_zmatrix_com(tmp_path):
    gaussian_input = tmp_path / "zmat.com"
    gaussian_input.write_text(
        """
#p hf/sto-3g geom=coord=zmat

old z-matrix input

0 1
C
H 1 R1

R1=1.09
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Cartesian coordinates"):
        read_geometry_input(gaussian_input)


def test_semiexperimental_reads_gaussian_cartesian_com_with_blank_title(tmp_path):
    gaussian_input = tmp_path / "blank_title.com"
    gaussian_input.write_text(
        """
#p hf/sto-3g


0 1
C 0.0 0.0 0.0
H 0.0 0.0 1.0
""",
        encoding="utf-8",
    )

    geometry = read_geometry_input(gaussian_input)

    assert geometry.comment == "blank_title"
    assert geometry.atoms == ("C", "H")


def test_semiexperimental_reads_standard_job_toml(tmp_path):
    observations = tmp_path / "isotopologues.toml"
    observations.write_text(
        """
[[isotopologues]]
label = "parent"
[isotopologues.constants]
A_MHz = 1000.0
B_MHz = 800.0
C_MHz = 600.0
""",
        encoding="utf-8",
    )
    job_path = tmp_path / "water.mse.toml"
    job_path.write_text(
        f"""
schema = "{SEMIEXP_JOB_SCHEMA}"
title = "water SE fit"

[files]
observations = "{observations.name}"

[fit]
observable = "moments"
rotational_components = "auto"
max_step = 0.2
prune_condition = 0.0

[geometry]
units = "angstrom"
atoms = [
  ["O", 0.000000, 0.000000, 0.000000],
  ["H", 0.000000, 0.000000, 0.957200],
  ["H", 0.926600, 0.000000, -0.239600],
]

[constraints]
fix_hydrogen_parameters = true
modredundant = [
  "B 1 2 F",
]
fixed_gic_patterns = ["angle(2,1,3)"]

[[qm_predicates]]
pattern = "GIC001"
value = 1.0
sigma = 0.1
source = "test"

[[parameter_classes]]
name = "OH"
mode = "shared"
patterns = ["bond(1,2)", "bond(1,3)"]
""",
        encoding="utf-8",
    )

    job = read_semiexperimental_job(job_path)
    geometry = read_geometry_input(job_path)

    assert job.observations == observations
    assert job.geometry.atoms == ("O", "H", "H")
    assert geometry.source_format == "merlino_semiexp_job"
    assert "bond(1,2)" in job.fixed_parameters
    assert "angle(2,1,3)" in job.fixed_parameters
    assert HYDROGEN_PARAMETER_CONSTRAINT in job.fixed_parameters
    assert job.qm_predicates[0].source == "test"
    assert job.parameter_classes[0].name == "OH"


def test_semiexperimental_rejects_dummy_atoms_in_cartesian_com(tmp_path):
    gaussian_input = tmp_path / "dummy.com"
    gaussian_input.write_text(
        """
#p hf/sto-3g

dummy atom

0 1
X 0.0 0.0 0.0
C 0.0 0.0 1.0
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Invalid atom token"):
        read_geometry_input(gaussian_input)


def test_semiexperimental_substitution_parser():
    assert parse_substitutions("2:13;5:18;6:D;7:T") == {2: 13, 5: 18, 6: 2, 7: 3}
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
    geometry_input = tmp_path / "water_initial.com"
    geometry_input.write_text(
        "\n".join(
            [
                "#p hf/sto-3g opt=modredundant",
                "",
                "distorted water",
                "",
                "0 1",
                *[f"{atom} {x:.8f} {y:.8f} {z:.8f}" for atom, (x, y, z) in zip(atoms, initial)],
                "",
            ]
        ),
        encoding="utf-8",
    )
    parent_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, target)))
    d1_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, target, [None, 2, None])))
    observations = (
        IsotopologueObservation("parent", parent_constants),
        IsotopologueObservation("D1", d1_constants, substitutions={2: 2}),
    )
    request = SemiexperimentalFitRequest(geometry_input, observations)
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
        "line_search_stalled",
        "max_iter",
    }
    assert all(np.isfinite(parameter.sigma) for parameter in result.parameters)
    assert any(item.kind == "bond" and item.value_angstrom is not None for item in result.geometry_parameters)
    assert any(item.kind == "bond" and item.sigma_angstrom is not None for item in result.geometry_parameters)
    assert any(item.kind == "angle" and item.value_degree is not None for item in result.geometry_parameters)
    assert any(item.kind == "angle" and item.sigma_degree is not None for item in result.geometry_parameters)
    assert len(result.rotational_constants) == 3 * len(observations)
    assert any(item.component == "A" for item in result.rotational_constants)
    assert (tmp_path / "semiexp" / "semiexp_geometry.xyz").exists()
    assert (tmp_path / "semiexp" / "semiexp_parameters.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_geometry_parameters.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_residuals.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_rotational_constants.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_covariance.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_correlation.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_hessian.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_hessian_eigenvalues.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_diagnostics.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_influence.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_high_correlations.csv").exists()
    influence_text = (tmp_path / "semiexp" / "semiexp_influence.csv").read_text(encoding="utf-8")
    assert "chi_square_contribution" in influence_text
    diagnostics_text = (tmp_path / "semiexp" / "semiexp_diagnostics.csv").read_text(encoding="utf-8")
    assert "incremental_rank" in diagnostics_text
    assert "parameter_scale_min" in diagnostics_text
    assert (tmp_path / "semiexp" / "semiexp_manifest.json").exists()
    rotconst_text = (tmp_path / "semiexp" / "semiexp_rotational_constants.csv").read_text(encoding="utf-8")
    assert "corrected_experimental_MHz" in rotconst_text
    assert "difference_MHz" in rotconst_text


def test_semiexperimental_fit_can_use_hessian_free_symmetry_cartesians(tmp_path):
    fchk = Path("gui/tests/gaussian/h2o.fchk")
    hessian_input = hessian_input_from_gaussian_fchk(fchk)
    atoms = ("H", "O", "H")
    coords = hessian_input.cartesian_coordinates_bohr * BOHR_TO_ANGSTROM
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water from fchk",
                *[f"{atom} {x:.12f} {y:.12f} {z:.12f}" for atom, (x, y, z) in zip(atoms, coords)],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords))),
    )
    model = cartesian_symmetry_coordinate_model(atoms, coords)

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(
            xyz,
            (observation,),
            coordinate_model="cartesian_symmetry",
        ),
        max_iter=1,
        outdir=tmp_path / "cartesian_symmetry",
    )
    manifest = json.loads((tmp_path / "cartesian_symmetry" / "semiexp_manifest.json").read_text(encoding="utf-8"))
    report = (tmp_path / "cartesian_symmetry" / "semiexp_report.txt").read_text(encoding="utf-8")

    assert model.point_group == "C2v"
    assert model.model_kind == "cartesian_symmetry"
    assert model.cartesian_from_q.shape == (3 * len(atoms), 3 * len(atoms) - 6)
    assert "A1" in set(model.irreps)
    assert result.diagnostics.coordinate_model == "cartesian_symmetry"
    assert result.b_matrix.shape == (len(result.gic_labels), 3 * len(atoms))
    assert all(("irrep=A1" in parameter.name) == parameter.active for parameter in result.parameters)
    assert manifest["backend"]["coordinate_model"] == "symmetry-cartesian"
    assert manifest["parameters"]["coordinate_generation"]["active_subspace"] == "totally symmetric symmetry-adapted Cartesian displacements only"
    assert "coordinate_basis = totally symmetric Hessian-free symmetry-adapted Cartesian displacements" in report
    assert "hessian =" not in report


def test_semiexperimental_topological_dihedral_errors_are_propagated():
    from merlino_semiexp.fit import _geometry_parameters

    atoms = ("C", "C", "C", "C")
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [2.5, 1.0, 0.0],
            [3.5, 1.0, 1.0],
        ],
        dtype=float,
    )
    fit_prims = (
        Primitive("bond", (0, 1)),
        Primitive("bond", (1, 2)),
        Primitive("bond", (2, 3)),
        Primitive("angle", (0, 1, 2)),
        Primitive("angle", (1, 2, 3)),
        Primitive("dihedral", (0, 1, 2, 3)),
    )

    rows = _geometry_parameters(
        atoms,
        coords,
        fit_prims=fit_prims,
        fit_u_matrix=np.eye(len(fit_prims)),
        active_mask=np.ones(len(fit_prims), dtype=bool),
        transform=np.eye(len(fit_prims)),
        covariance=np.eye(len(fit_prims)) * 1.0e-6,
    )

    assert any(item.kind == "dihedral" and item.value_degree is not None for item in rows)
    assert any(item.kind == "dihedral" and item.sigma_degree is not None for item in rows)


def test_semiexperimental_weak_parameter_pruning_is_deterministic():
    from merlino_semiexp.fit import _weak_parameter_patterns

    weighted_jac = np.diag([10.0, 1.0, 0.01])
    names = (
        "GIC001 GICForge A1Str0001 irrep=A1",
        "GIC002 GICForge A1Ang0001 irrep=A1",
        "GIC003 GICForge A1Ang0002 irrep=A1",
    )

    assert _weak_parameter_patterns(names, weighted_jac, 20.0) == ("A1Ang0002",)


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


def test_semiexperimental_primitive_constraints_do_not_disable_containing_gics(tmp_path):
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
        SemiexperimentalFitRequest(xyz, (observation,), fixed_parameters=("bond(1,2)",)),
        max_iter=1,
    )

    active_count = sum(parameter.active for parameter in result.parameters)
    assert any("bond(1,2)" in parameter.name and parameter.active for parameter in result.parameters)
    assert result.jacobian.shape[1] < active_count


def test_semiexperimental_primitive_constraints_expand_by_symmetry(tmp_path):
    atoms = ("O", "H", "H")
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    prims, _u_matrix, _labels = _gic_model(
        coords,
        z_numbers,
        backend=_make_gicforge_backend(atoms, tmp_path),
    )
    fixed = _fixed_primitives_from_patterns(("bond(1,2)",))

    expanded = _symmetry_expanded_fixed_primitives(atoms, coords, prims, fixed)
    expanded_keys = {_primitive_constraint_key(primitive) for primitive in expanded}

    assert expanded_keys == {
        ("bond", (0, 1), 0),
        ("bond", (0, 2), 0),
    }


def test_semiexperimental_hydrogen_constraint_generates_h_primitives(tmp_path):
    atoms = ("O", "H", "H")
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    prims, _u_matrix, _labels = _gic_model(
        coords,
        z_numbers,
        backend=_make_gicforge_backend(atoms, tmp_path),
    )

    fixed = _hydrogen_fixed_primitives(atoms, prims, (HYDROGEN_PARAMETER_CONSTRAINT,))
    fixed_keys = {_primitive_constraint_key(primitive) for primitive in fixed}

    assert fixed
    assert len(fixed_keys) == 3
    assert all(any(atoms[atom] == "H" for atom in primitive.atoms) for primitive in fixed)
    assert ("bond", (0, 1), 0) in fixed_keys
    assert ("bond", (0, 2), 0) in fixed_keys
    assert ("angle", (1, 0, 2), 0) in fixed_keys


def test_semiexperimental_parameter_classes_share_and_fix_parameters(tmp_path):
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
        SemiexperimentalFitRequest(
            xyz,
            (observation,),
            parameter_classes=(
                ParameterClassConstraint("OH_stretches", ("bond(1,2)", "bond(1,3)"), "shared"),
                ParameterClassConstraint("antisymmetric_stretch", ("GIC002",), "fixed"),
            ),
        ),
        max_iter=1,
        outdir=tmp_path / "semiexp_classes",
    )

    shared = [parameter for parameter in result.parameters if parameter.parameter_class == "OH_stretches"]
    fixed = [parameter for parameter in result.parameters if parameter.parameter_class == "antisymmetric_stretch"]

    assert shared
    assert all(parameter.active for parameter in shared)
    assert fixed and all(not parameter.active for parameter in fixed)
    assert result.jacobian.shape[1] == 1
    params_text = (tmp_path / "semiexp_classes" / "semiexp_parameters.csv").read_text(encoding="utf-8")
    assert "parameter_class" in params_text
    assert "OH_stretches" in params_text


def test_semiexperimental_kraitchman_comparison_for_single_substitution(tmp_path):
    atoms = ["O", "H", "H"]
    coords = np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 0.9572],
            [0.9266, 0.0000, -0.2396],
        ],
        dtype=float,
    )
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(["3", "water", *[f"{atom} {x:.8f} {y:.8f} {z:.8f}" for atom, (x, y, z) in zip(atoms, coords)]])
        + "\n",
        encoding="utf-8",
    )
    observations = (
        IsotopologueObservation("parent", RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords)))),
        IsotopologueObservation(
            "D1",
            RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords, [None, 2, None]))),
            substitutions={2: 2},
        ),
    )
    direct_rows = kraitchman_comparison(atoms, coords, observations)
    direct_seed = kraitchman_seed_geometry(atoms, coords, observations, direct_rows)

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(xyz, observations),
        max_iter=1,
        outdir=tmp_path / "semiexp_krai",
    )

    assert direct_seed is not None
    assert direct_seed.method == "direct_substitution"
    assert direct_seed.fitted_atom_indices == (2,)
    assert direct_seed.coordinates_angstrom[1, 0] == pytest.approx(direct_rows[0].signed_kraitchman_angstrom)
    assert len(result.kraitchman) == 3
    assert result.kraitchman_seed is not None
    assert {row.coordinate for row in result.kraitchman} == {"a", "b", "c"}
    assert all(row.isotopologue == "D1" for row in result.kraitchman)
    assert (tmp_path / "semiexp_krai" / "semiexp_kraitchman.csv").exists()
    assert (tmp_path / "semiexp_krai" / "semiexp_kraitchman_geometry.xyz").exists()


def test_semiexperimental_gic_preview_and_html_report(tmp_path):
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9572], [0.9266, 0.0, -0.2396]])
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(["3", "water", *[f"{a} {x:.8f} {y:.8f} {z:.8f}" for a, (x, y, z) in zip(atoms, coords)]])
        + "\n",
        encoding="utf-8",
    )
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords))),
    )
    request = SemiexperimentalFitRequest(xyz, (observation,))

    preview = preview_semiexperimental_gics(xyz, request.observations)
    result = fit_semiexperimental_geometry(request, max_iter=1)
    report = write_semiexperimental_html_report(tmp_path / "semiexp_report.html", result, request)

    assert "Non-redundant GICs" in preview.text
    assert preview.gic_labels
    assert preview.rows
    assert {row.kind for row in preview.rows}.issubset({"bond", "angle", "dihedral", "out_of_plane", "linear_bend", "ring", "mixed"})
    assert any(item.mode in {"shared", "fixed"} for item in preview.suggested_classes)
    report_text = report.read_text(encoding="utf-8")
    assert "Merlino Semiexperimental Geometry Report" in report_text
    assert "Final Cartesian Geometry Parameters" in report_text
    assert "Rotational Constants" in report_text
    assert "Corrected experimental / MHz" in report_text
    assert "Angle or dihedral / degree" in report_text
    tables = semiexperimental_latex_tables(result)
    assert {"parameters", "rotational_constants", "residuals", "kraitchman"} == set(tables)
    assert "\\begin{tabular}" in tables["parameters"]


def test_semiexperimental_validation_and_conditioning_preview(tmp_path):
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9572], [0.9266, 0.0, -0.2396]])
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(["3", "water", *[f"{a} {x:.8f} {y:.8f} {z:.8f}" for a, (x, y, z) in zip(atoms, coords)]])
        + "\n",
        encoding="utf-8",
    )
    observations = (
        IsotopologueObservation("parent", RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords)))),
        IsotopologueObservation(
            "D1",
            RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords, [None, 2, None]))),
            substitutions={2: 2},
        ),
    )
    request = SemiexperimentalFitRequest(
        xyz,
        observations,
        parameter_classes=(ParameterClassConstraint("OH_stretches", ("bond(1,2)", "bond(1,3)"), "shared"),),
    )

    issues = validate_semiexperimental_request(request)
    conditioning = preview_semiexperimental_conditioning(request)

    assert not [item for item in issues if item.severity == "error"]
    assert conditioning.n_observations >= 3
    assert conditioning.n_effective_parameters >= 1
    assert "condition number" in conditioning.text


def test_semiexperimental_gic_preview_keeps_angstrom_topology_for_cyclopentadiene(tmp_path):
    xyz = tmp_path / "cyclopentadiene.xyz"
    xyz.write_text(
        "\n".join(
            [
                "11",
                "cyclopentadiene MSR geometry",
                "C 0.0000000648 1.21099405 -0.0000437769",
                "C 1.16966652 0.34210493 0.0000360217",
                "C -1.16966649 0.34210505 0.0000319263",
                "C 0.731740874 -0.975616979 -0.0000177892",
                "C -0.731740976 -0.975616904 -0.0000203515",
                "H 2.26478781 0.374863776 0.0001274146",
                "H -2.26478777 0.374864010 0.0001200854",
                "H 1.26945754 -1.92001108 0.0000062677",
                "H -1.26945774 -1.92001095 0.0000018591",
                "H 0.0000585394 1.87871494 -0.8689120000",
                "H -0.0000583448 1.87871720 0.8688227060",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    preview = preview_semiexperimental_gics(xyz)

    assert preview.gic_labels
    assert any("dihedral" in label for label in preview.gic_labels)


def test_semiexperimental_fit_uses_adaptive_gicforge_model(tmp_path, monkeypatch):
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9572], [0.9266, 0.0, -0.2396]])
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(["3", "", *[f"{a} {x:.8f} {y:.8f} {z:.8f}" for a, (x, y, z) in zip(atoms, coords)]])
        + "\n",
        encoding="utf-8",
    )
    calls = []

    def fake_run_gicforge(workdir):
        workdir = Path(workdir)
        xyzin_lines = (workdir / "xyzin").read_text(encoding="utf-8").splitlines()
        assert xyzin_lines[1] == ""
        (workdir / "provout").write_text(" Point Group from symm.f: C2v\n", encoding="utf-8")
        gauin = workdir / "gauin.symm"
        lines = [
            " A1Str0001=[ 0.7071*R(  1,  2)+0.7071*R(  1,  3)]",
            " B2Str0001=[ 0.7071*R(  1,  2)-0.7071*R(  1,  3)]",
            " A1Ang0001=[ 1.00000*A(  2,  1,  3)]",
            " B1Lin0001 = L(  2,  1,  3,  0, -1)",
            " A2Oop0001 = U(  2,  1,  3,  2)",
        ]
        gauin.write_text(
            "\n".join(
                lines
            )
            + "\n",
            encoding="utf-8",
        )
        (workdir / "gicsym").write_text(
            "\n".join(
                [
                    "name,irrep",
                    "A1Str0001,A1",
                    "B2Str0001,B2",
                    "A1Ang0001,A1",
                    "B1Lin0001,B1",
                    "A2Oop0001,A2",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        calls.append(workdir)
        return SimpleNamespace(files={"gauin.symm": gauin, "gicsym": workdir / "gicsym"})

    monkeypatch.setattr("merlino_semiexp.fit.run_gicforge", fake_run_gicforge)
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords))),
    )

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(xyz, (observation,)),
        max_iter=1,
        outdir=tmp_path / "run",
    )

    manifest = json.loads((tmp_path / "run" / "semiexp_manifest.json").read_text(encoding="utf-8"))
    assert calls
    assert len(calls) == 1
    assert all("GICForge" in parameter.name for parameter in result.parameters)
    assert any(parameter.active for parameter in result.parameters)
    assert any(not parameter.active for parameter in result.parameters)
    assert any("B1Lin" in label for label in result.gic_labels)
    assert any("A2Oop" in label for label in result.gic_labels)
    assert result.b_matrix.shape[0] == 5
    assert manifest["backend"]["coordinate_model"] == "gicforge-frozen-definition"
    assert "run once" in manifest["parameters"]["coordinate_generation"]["primitive_source"]
    assert "frozen GIC schema" in manifest["parameters"]["coordinate_generation"]["line_search"]
    assert "secant-updated B projector" in manifest["parameters"]["coordinate_generation"]["line_search"]
    assert manifest["parameters"]["gicforge_calls"] == 1
    assert manifest["parameters"]["b_projector_analytic_refreshes"] >= 1
    assert manifest["parameters"]["coordinate_generation"]["active_subspace"] == "GICForge-assigned A1 coordinates only"


def test_semiexperimental_validation_flags_bad_classes_and_isotopes(tmp_path):
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
    obs = IsotopologueObservation(
        "bad",
        RotationalConstants(1000.0, 800.0, 600.0),
        substitutions={1: 2},
    )
    request = SemiexperimentalFitRequest(
        xyz,
        (obs,),
        parameter_classes=(ParameterClassConstraint("missing", ("not_a_gic",), "shared"),),
    )

    issues = validate_semiexperimental_request(request)

    assert any("deuterium substitution on non-H" in item.message for item in issues)
    assert any("matches no GIC" in item.message for item in issues)


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
