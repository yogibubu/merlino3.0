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
from merlino_semiexp import (
    CorrectedRotationalConstants,
    DEFAULT_SEMIEXP_OBSERVABLE,
    DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
    ElectronicCorrection,
    IsotopologueObservation,
    ParameterClassConstraint,
    QMParameterPredicate,
    RotationalConstants,
    SemiexperimentalFitRequest,
    VibrationalCorrection,
    corrected_constants_rows,
    fit_semiexperimental_geometry,
    parse_substitutions,
    preview_semiexperimental_conditioning,
    preview_semiexperimental_gics,
    read_observations,
    read_observations_csv,
    semiexperimental_latex_tables,
    validate_semiexperimental_request,
    write_semiexperimental_html_report,
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
    electronic = ElectronicCorrection(0.25, 0.5, -0.25, source="electronic")
    corrected = CorrectedRotationalConstants(observed, correction, electronic).equilibrium

    assert corrected.as_tuple() == (998.75, 801.5, 599.75)


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
    assert not any(parameter.active for parameter in result.parameters)
    assert result.diagnostics.convergence_reason == "no_active_totally_symmetric_parameters"


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

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(xyz, observations),
        max_iter=1,
        outdir=tmp_path / "semiexp_krai",
    )

    assert len(result.kraitchman) == 3
    assert {row.coordinate for row in result.kraitchman} == {"a", "b", "c"}
    assert all(row.isotopologue == "D1" for row in result.kraitchman)
    assert (tmp_path / "semiexp_krai" / "semiexp_kraitchman.csv").exists()


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
    assert "Merlino Semiexperimental Geometry Report" in report.read_text(encoding="utf-8")
    tables = semiexperimental_latex_tables(result)
    assert {"parameters", "residuals", "kraitchman"} == set(tables)
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


def test_semiexperimental_fit_always_uses_iterative_gicforge(tmp_path, monkeypatch):
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
        gauin = workdir / "gauin"
        gauin.write_text(
            "\n".join(
                [
                    " Stre0001=[ 0.7071*R(  1,  2)+0.7071*R(  1,  3)]",
                    " Stre0002=[ 0.7071*R(  1,  2)-0.7071*R(  1,  3)]",
                    " SymD0001 =[ 1.00000*A(  2,  1,  3)]",
                    " LAng0001 = L(  2,  1,  3,  0, -1)",
                    " OuPl0001 = U(  2,  1,  3,  2)",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        calls.append(workdir)
        return SimpleNamespace(files={"gauin": gauin})

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
    assert all("GICForge" in parameter.name for parameter in result.parameters)
    assert any(parameter.active for parameter in result.parameters)
    assert any(not parameter.active for parameter in result.parameters)
    assert any("LAng" in label for label in result.gic_labels)
    assert any("OuPl" in label for label in result.gic_labels)
    assert result.b_matrix.shape[0] == 5
    assert manifest["backend"]["coordinate_model"] == "gicforge-iterative-readallgic"
    assert "GICForge ReadAllGIC" in manifest["parameters"]["coordinate_generation"]["primitive_source"]
    assert manifest["parameters"]["coordinate_generation"]["active_subspace"] == "totally symmetric GICForge coordinates only"


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
