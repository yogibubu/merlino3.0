from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np

from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_core import ScientificValidationError
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
    finite_difference_constraint_b_matrix,
    mass_vector_for_observation,
    parse_gaussian_style_constraints,
    corrected_constants_rows,
    fit_semiexperimental_geometry,
    kraitchman_comparison,
    kraitchman_seed_geometry,
    parse_substitutions,
    preview_semiexperimental_conditioning,
    preview_semiexperimental_gics,
    read_geometry_input,
    read_msr_legacy_input,
    read_observations,
    read_observations_csv,
    read_semiexperimental_job,
    semiexperimental_latex_tables,
    validate_semiexperimental_request,
    write_semiexperimental_html_report,
    write_observations_csv,
    is_msr_legacy_file,
)
from merlino_vpt2_vci.gaussian_qff import hessian_input_from_gaussian_fchk
from merlino_core import repo_root
from merlino_semiexp.fit import (
    MeasurementModel,
    SemiexperimentalFitDiagnostics,
    SemiexperimentalGeometryParameter,
    SemiexperimentalParameter,
    _combined_primitive_constraint_b_matrix,
    _atomic_number,
    _auto_resolve_isotopic_substitutions,
    _dynamic_parameter_scales,
    _finite_difference_constraint_b_matrix,
    _fixed_primitives_from_patterns,
    _gic_expression_definitions_from_patterns,
    _gic_fixed_patterns,
    _gic_expression_constraint_targets,
    _gic_expression_constraint_values,
    _gic_expression_constraints_from_patterns,
    _gic_model,
    _hydrogen_fixed_primitives,
    _isotopic_mapping_warning_rows,
    _linear_primitive_constraints_from_patterns,
    _make_gicforge_backend,
    _primitive_constraint_key,
    _rank_revealing_lm_step,
    _rotational_residual_stats,
    _robust_sqrt_weights,
    _semiexp_warning_rows,
    _stationary_point_type,
    _svd_diagnostics_csv,
    _symmetry_expanded_fixed_primitives,
    _uncertainty_diagnostics_csv,
    _warnings_csv,
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

    robust = SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs,), robust_loss="huber", robust_scale=2.0)
    robust.validate()

    with pytest.raises(ValueError):
        SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs,), robust_loss="bad").validate()

    with pytest.raises(ValueError):
        SemiexperimentalFitRequest(tmp_path / "geom.xyz", (obs,), robust_scale=-1.0).validate()


def test_semiexperimental_robust_weights_apply_only_to_experimental_rows():
    residual = np.array([0.1, 20.0, 0.2, 50.0], dtype=float)
    weights, scale, downweighted_rows, downweighted_isotopologues = _robust_sqrt_weights(
        residual,
        "huber",
        1.0,
        experimental_rows=3,
        row_groups=((0, 1), (2,)),
    )

    assert scale == pytest.approx(1.0)
    assert downweighted_rows == 2
    assert downweighted_isotopologues == 1
    assert weights[0] == pytest.approx(weights[1])
    assert weights[0] < 1.0
    assert weights[1] < 1.0
    assert weights[2] == pytest.approx(1.0)
    assert weights[3] == pytest.approx(1.0)


def test_semiexperimental_dynamic_column_scaling_equilibrates_jacobian():
    jac = np.array([[1.0e3, 1.0], [2.0e3, 2.0], [3.0e3, 3.0]], dtype=float)
    scales = _dynamic_parameter_scales(jac, np.ones(2, dtype=float))
    before = np.linalg.norm(jac, axis=0)
    after = np.linalg.norm(jac * scales[None, :], axis=0)

    assert scales[0] < scales[1]
    assert before[0] / before[1] > 1.0e2
    assert after[0] / after[1] == pytest.approx(1.0)


def test_semiexperimental_rank_revealing_step_handles_dependent_columns():
    jac = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]], dtype=float)
    residual = np.array([1.0, 2.0, 3.0], dtype=float)
    step = _rank_revealing_lm_step(jac, residual, damping=1.0e-8)
    predicted = residual - jac @ step

    assert np.all(np.isfinite(step))
    assert step[0] == pytest.approx(step[1])
    assert float(predicted @ predicted) < float(residual @ residual)


def test_semiexperimental_svd_diagnostics_report_near_null_combinations():
    labels = ("q1", "q2")
    jac = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]], dtype=float)
    text = _svd_diagnostics_csv(labels, jac)

    assert "dominant_coordinate_combination" in text
    assert "q1" in text
    assert "q2" in text
    assert ",1," in text


def test_semiexperimental_diagnostic_warnings_are_machine_readable():
    diagnostics = SemiexperimentalFitDiagnostics(
        convergence_reason="max_iter",
        objective=1.0,
        weighted_rms=1.0,
        reduced_chi_square=1.0,
        rank=1,
        incremental_rank=1,
        condition_number=1.0e12,
        damping=1.0e-6,
        accepted_steps=1,
        rejected_steps=0,
        max_iterations=1,
        n_optimized_parameters=2,
        observable="moments",
        components=("Ia", "Ib"),
        planar=True,
        robust_loss="cauchy",
        robust_scale=1.0,
        robust_downweighted_observations=2,
        robust_downweighted_isotopologues=1,
    )
    model = MeasurementModel(
        observable="moments",
        components=("Ia", "Ib"),
        labels=(("iso_low", "Ia"), ("iso_low", "Ib")),
        observed=np.array([1.0, 2.0]),
        weights=np.ones(2),
        n_experimental_rows=2,
        planar=True,
    )
    parameters = (
        SemiexperimentalParameter("q1", 0.0, 1.0e-3, True),
        SemiexperimentalParameter("q2", 0.0, 1.0, True),
    )
    geometry = (
        SemiexperimentalGeometryParameter(
            "bond",
            "R(1,2)",
            (1, 2),
            ("C", "C"),
            value_angstrom=1.40,
            sigma_angstrom=0.02,
        ),
    )
    jac = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]], dtype=float)

    rows = _semiexp_warning_rows(
        diagnostics,
        ("q1", "q2"),
        parameters,
        geometry,
        jac,
        model,
        np.array([0.2, 0.2]),
    )
    codes = {row.code for row in rows}
    csv_text = _warnings_csv(rows)

    assert "rank_deficient" in codes
    assert "small_singular_value" in codes
    assert "planar_pair_ill_conditioned" in codes
    assert "low_robust_isotopologue_weight" in codes
    assert "large_geometry_uncertainty" in codes
    assert "severity,code,message,context" in csv_text
    assert "iso_low" in csv_text


def test_semiexperimental_flags_suspicious_isotopologue_atom_mapping():
    atoms = ("O", "C", "O", "H")
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.25, 0.0, 0.0],
            [2.55, 0.45, 0.0],
            [2.95, 1.15, 0.35],
        ],
        dtype=float,
    )
    parent_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords)))
    right_o18_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords, [None, None, 18, None])))
    observations = (
        IsotopologueObservation("parent", parent_constants),
        IsotopologueObservation("O18_wrong_atom", right_o18_constants, substitutions={1: 18}),
    )

    rows = _isotopic_mapping_warning_rows(atoms, coords, observations)

    assert any(row.code == "isotopologue_mapping_suspicious" for row in rows)
    assert any("input_atom=1" in row.context and "suggested_atom=3" in row.context for row in rows)


def test_semiexperimental_autocorrects_clear_single_isotopologue_mapping():
    atoms = ("O", "C", "O", "H")
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.25, 0.0, 0.0],
            [2.55, 0.45, 0.0],
            [2.95, 1.15, 0.35],
        ],
        dtype=float,
    )
    parent_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords)))
    right_o18_constants = RotationalConstants(*rotational_constants_MHz(_structure(atoms, coords, [None, None, 18, None])))
    observations = (
        IsotopologueObservation("parent", parent_constants),
        IsotopologueObservation("O18_wrong_atom", right_o18_constants, substitutions={1: 18}),
    )

    resolved, warnings = _auto_resolve_isotopic_substitutions(atoms, coords, observations)

    assert resolved[1].substitutions == {3: 18}
    assert any(row.code == "isotopologue_mapping_autocorrected" for row in warnings)
    assert any("input_atom=1" in row.context and "used_atom=3" in row.context for row in warnings)


def test_semiexperimental_stationary_point_uses_numerical_eigenvalue_tolerance():
    assert _stationary_point_type(np.array([5.8e4, 1.0e8, 6.9e12])) == "minimum"
    assert _stationary_point_type(np.array([-1.0e-2, 1.0, 10.0])) == "transition_state_or_saddle"


def test_semiexperimental_rotational_residual_stats_include_mean_square():
    rows = (
        SimpleNamespace(difference_MHz=0.1),
        SimpleNamespace(difference_MHz=-0.3),
    )

    nrows, rms, mean_square, scaled_mean_square, max_abs = _rotational_residual_stats(rows)

    assert nrows == 2
    assert rms == pytest.approx((0.05) ** 0.5)
    assert mean_square == pytest.approx(0.05)
    assert scaled_mean_square == pytest.approx(50.0)
    assert max_abs == pytest.approx(0.3)


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
QFIX=[R(1,3)-R(1,2)] Value=0.0
HOH(Frozen)=A(2,1,3)
""",
        encoding="utf-8",
    )

    geometry = read_geometry_input(gaussian_input)

    assert geometry.source_format == "gaussian_cartesian_modredundant"
    assert geometry.comment == "water constrained"
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.coordinates_angstrom.shape == (3, 3)
    assert "R(1,2) Frozen" in geometry.fixed_parameters
    assert "A(2,1,3) Frozen" in geometry.fixed_parameters
    assert "QFIX=[R(1,3)-R(1,2)] Value=0.0" in geometry.fixed_parameters
    assert "HOH(Frozen)=A(2,1,3)" in geometry.fixed_parameters


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
fixed_gic_patterns = ["A(2,1,3) Frozen"]
gic_constraints = ["QFIX=[GIC001+2*GIC002] Value=0.0"]

[[qm_predicates]]
pattern = "GIC001"
value = 1.0
sigma = 0.1
source = "test"

[[parameter_classes]]
name = "OH"
mode = "shared"
patterns = ["R(1,2)", "R(1,3)"]
""",
        encoding="utf-8",
    )

    job = read_semiexperimental_job(job_path)
    geometry = read_geometry_input(job_path)

    assert job.observations == observations
    assert job.geometry.atoms == ("O", "H", "H")
    assert geometry.source_format == "merlino_semiexp_job"
    assert "R(1,2) Frozen" in job.fixed_parameters
    assert "A(2,1,3) Frozen" in job.fixed_parameters
    assert "QFIX=[GIC001+2*GIC002] Value=0.0" in job.fixed_parameters
    assert HYDROGEN_PARAMETER_CONSTRAINT in job.fixed_parameters
    assert job.qm_predicates[0].source == "test"
    assert job.parameter_classes[0].name == "OH"


def test_semiexperimental_reads_legacy_msr_input_without_required_labels(tmp_path):
    msr_path = tmp_path / "minimal.msr.inp"
    msr_path.write_text(
        """
#m optim=(method=gaun,coord=zmat) geom=coord=zmat
#d niso=2

C
X  1 #XX
H  1 #CH  2 #A90

CH = 1.0900
XX = 1.0000
A90 = 90.0

12.000000
 1.000000

12.000000
 2.000000

bexp
 1000.0 800.0 600.0
  900.0 700.0 500.0

dbvib
 1.0 2.0 3.0
 4.0 5.0 6.0

weights
 1.0 4.0 9.0
 16.0 25.0 36.0
""",
        encoding="utf-8",
    )

    legacy = read_msr_legacy_input(msr_path)
    geometry = read_geometry_input(msr_path)
    observations = read_observations(msr_path)

    assert is_msr_legacy_file(msr_path)
    assert legacy.geometry.source_format == "msr_legacy_zmatrix"
    assert geometry.atoms == ("C", "H")
    assert geometry.coordinates_angstrom.shape == (2, 3)
    assert "R(1,2) Frozen" in geometry.fixed_parameters
    assert observations[0].label == "parent"
    assert observations[1].label == "iso_002"
    assert observations[1].substitutions == {2: 2}
    assert observations[0].corrected.as_tuple() == pytest.approx((1001.0, 802.0, 603.0))
    assert observations[1].weights is not None
    assert observations[1].weights.as_tuple() == pytest.approx((16.0, 25.0, 36.0))


def test_legacy_msr_nitrobenzene_zmatrix_keeps_closed_ring_and_nonredundant_gics(tmp_path):
    msr_path = tmp_path / "nitrobenzene.msr.inp"
    msr_path.write_text(
        """
#m optim=(method=gaun,coord=zmat) geom=coord=zmat
#d niso=1

 C
 N  1  CN
 X  2  #XX  1  #A90
 X  3  #XX  2  #A90 1 #D180
 C  1  CC   2  CCN  3 #D000
 H  5  CH   1  CCH  2 #D000
 C  1  CC   2  CCN  3 #D180
 H  7  CH   1  CCH  2 #D000
 C  5  CC1  1  CCC  2 #D180
 H  9  CH1  5  CCH1 1 #D180
 C  7  CC1  1  CCC  2 #D180
 H 11  CH1  7  CCH1 1 #D180
 C  2  CC2  3  #A90 4 #D180
 H  2  CH2  3  #A90 4 #D180
 O  2  NO   1  CNO  5 #D000
 O  2  NO   1  CNO  5 #D180

  CN   =    1.4734
  NO   =    1.2245
  CC   =    1.389
  CH   =    1.0807
  CC1  =    1.3911
  CH1  =    1.083
  CC2  =    4.2267
  CH2  =    5.3101
  CNO  =   117.5927
  CCN  =   118.7407
  CCH  =   119.7612
  CCC  =   118.3592
 CCH1  =   119.5806
  XX   =    1.0
  A90  =    90.0
 D180  =   180.0
 D000  =     0.0

 12.000000  \\parent
 14.000000
 12.000000
  1.000000
 12.000000
  1.000000
 12.000000
  1.000000
 12.000000
  1.000000
 12.000000
  1.000000
 16.000000
 16.000000

bexp
   3968.0780   1286.9203   972.6605 \\parent
""",
        encoding="utf-8",
    )

    from merlino_fit.survibfit.primitives import build_primitives
    from merlino_fit.survibfit.transforms import build_u, _vibrational_projector_local, _vibrational_rank
    from merlino_fit.topology.pipeline import build_topology_objects
    from topology.elements import atomic_number

    geometry = read_geometry_input(msr_path)
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    z_numbers = [atomic_number(atom) for atom in geometry.atoms]
    _continuous, graph, ringset, _synthons, _aromaticity = build_topology_objects(coords, z_numbers)

    assert sorted((i + 1, j + 1) for i, j in graph.bonds) == [
        (1, 2),
        (1, 3),
        (1, 5),
        (2, 13),
        (2, 14),
        (3, 4),
        (3, 7),
        (5, 6),
        (5, 9),
        (7, 8),
        (7, 11),
        (9, 10),
        (9, 11),
        (11, 12),
    ]
    assert [atom + 1 for atom in ringset.rings[0].atoms] == [1, 3, 7, 11, 9, 5]

    primitives = build_primitives(graph, coords)
    u_matrix = build_u(primitives, coords, Z=z_numbers, ringset=ringset)
    projector = _vibrational_projector_local(coords)
    rank = np.linalg.matrix_rank(u_matrix.T @ b_matrix_analytic(primitives, coords) @ projector, tol=1.0e-7)

    assert u_matrix.shape[1] == _vibrational_rank(coords) == 36
    assert rank == 36


def test_semiexperimental_reads_legacy_msr_cartesian_geometry(tmp_path):
    msr_path = tmp_path / "minimal.msr"
    msr_path.write_text(
        """
#m symmetry=isotopes

C  0.000000  0.000000  0.000000
H  0.000000  0.000000  1.090000

12.000000  \\parent
 1.000000

12.000000
 2.000000

bexp
 1000.0 800.0 600.0  \\parent
  900.0 700.0 500.0

dbvib
 0.0 0.0 0.0
 0.0 0.0 0.0
""",
        encoding="utf-8",
    )

    geometry = read_geometry_input(msr_path)
    observations = read_observations(msr_path)

    assert geometry.source_format == "msr_legacy_cartesian"
    assert geometry.atoms == ("C", "H")
    assert geometry.fixed_parameters == ()
    assert observations[0].label == "parent"
    assert observations[1].label == "iso_002"
    assert observations[1].substitutions == {2: 2}


def test_semiexperimental_reads_legacy_msr_cartesian_modredundant_constraints(tmp_path):
    msr_path = tmp_path / "cartesian_constraints.msr"
    msr_path.write_text(
        """
#m symmetry=isotopes

O  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.957200
H  0.926600  0.000000 -0.239600

constraints
B 1 2 F
A 2 1 3 F
end

16.000000  \\parent
 1.000000
 1.000000

16.000000
 2.000000
 1.000000

bexp
 1000.0 800.0 600.0  \\parent
  900.0 700.0 500.0
""",
        encoding="utf-8",
    )

    geometry = read_geometry_input(msr_path)
    observations = read_observations(msr_path)

    assert geometry.source_format == "msr_legacy_cartesian"
    assert geometry.atoms == ("O", "H", "H")
    assert "R(1,2) Frozen" in geometry.fixed_parameters
    assert "A(2,1,3) Frozen" in geometry.fixed_parameters
    assert observations[1].substitutions == {2: 2}


def test_semiexperimental_reads_robust_legacy_msr_zmatrix_variants(tmp_path):
    msr_path = tmp_path / "robust_zmatrix.msr.inp"
    msr_path.write_text(
        """
#m optim=(method=gaun,coord=zmat) geom=coord=zmat
symmetry=isotopes

C1
H2, 1, #RCH
H3  1  RCH  2  AHH

RCH  1.090000D+00
AHH  =  109.500000D+00

12.000000  \\parent
 1.000000
 1.000000

12.000000
 2.000000
 1.000000

bexp
 1000.0 800.0 600.0  \\parent
  900.0 700.0 500.0
""",
        encoding="utf-8",
    )

    geometry = read_geometry_input(msr_path)

    assert geometry.source_format == "msr_legacy_zmatrix"
    assert geometry.atoms == ("C", "H", "H")
    assert np.isfinite(geometry.coordinates_angstrom).all()
    assert geometry.coordinates_angstrom.shape == (3, 3)
    assert "R(1,2) Frozen" in geometry.fixed_parameters
    assert "R(1,3) Frozen" not in geometry.fixed_parameters


def test_semiexperimental_rejects_forward_zmatrix_references(tmp_path):
    msr_path = tmp_path / "bad_ref.msr.inp"
    msr_path.write_text(
        """
#m geom=coord=zmat

C
H 2 RCH

RCH = 1.09

12.000000
1.000000

bexp
1000.0 800.0 600.0
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="previous atoms"):
        read_geometry_input(msr_path)


def test_semiexperimental_does_not_treat_generic_inp_as_msr(tmp_path):
    generic = tmp_path / "minimal.inp"
    generic.write_text("C\nH 1 R\nR = 1.0\n", encoding="utf-8")

    assert not is_msr_legacy_file(generic)
    with pytest.raises(ValueError, match="Semiexperimental geometry input"):
        read_geometry_input(generic)
    with pytest.raises(ValueError, match="Semiexp observations"):
        read_observations(generic)


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
    request = SemiexperimentalFitRequest(geometry_input, observations, leave_one_out=True)
    initial_rms = _rotconst_rms(atoms, initial, observations)

    result = fit_semiexperimental_geometry(request, max_iter=8, outdir=tmp_path / "semiexp")

    assert result.rms_MHz < initial_rms
    assert result.diagnostics.observable == "moments"
    assert result.diagnostics.planar is True
    assert len(result.diagnostics.components) == 2
    assert set(result.diagnostics.components).issubset({"Ia", "Ib", "Ic"})
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
    assert len(result.leave_one_out) == len(observations)
    assert all(row.training_isotopologues == 1 for row in result.leave_one_out)
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
    assert (tmp_path / "semiexp" / "semiexp_svd_diagnostics.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_uncertainty_diagnostics.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_iteration_trace.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_constraints.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_warnings.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_leave_one_out.csv").exists()
    assert (tmp_path / "semiexp" / "semiexp_checkpoint.json").exists()
    influence_text = (tmp_path / "semiexp" / "semiexp_influence.csv").read_text(encoding="utf-8")
    assert "chi_square_contribution" in influence_text
    diagnostics_text = (tmp_path / "semiexp" / "semiexp_diagnostics.csv").read_text(encoding="utf-8")
    assert "incremental_rank" in diagnostics_text
    assert "parameter_scale_min" in diagnostics_text
    assert "robust_downweighted_isotopologues" in diagnostics_text
    assert (tmp_path / "semiexp" / "semiexp_manifest.json").exists()
    rotconst_text = (tmp_path / "semiexp" / "semiexp_rotational_constants.csv").read_text(encoding="utf-8")
    assert "corrected_experimental_MHz" in rotconst_text
    assert "difference_MHz" in rotconst_text
    report_text = (tmp_path / "semiexp" / "semiexp_report.txt").read_text(encoding="utf-8")
    assert "rotational_mean_square_MHz2" in report_text
    assert "[iteration_trace]" in report_text
    trace_text = (tmp_path / "semiexp" / "semiexp_iteration_trace.csv").read_text(encoding="utf-8")
    assert "objective_before" in trace_text
    assert "constraint_max_abs" in trace_text
    assert len(result.iteration_trace) >= 1


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
    assert result.diagnostics.components == ("A", "B")


def test_planar_moments_use_two_independent_components_and_explicit_pairs(tmp_path):
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

    auto = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(xyz, (observation,), rotational_components="auto"),
        max_iter=1,
    )
    explicit = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(xyz, (observation,), rotational_components="AC"),
        max_iter=1,
    )

    assert auto.diagnostics.planar is True
    assert auto.diagnostics.components == ("Ia", "Ib")
    assert explicit.diagnostics.components == ("Ia", "Ic")

    with pytest.raises(ScientificValidationError, match="ABC is redundant"):
        fit_semiexperimental_geometry(
            SemiexperimentalFitRequest(xyz, (observation,), rotational_components="ABC"),
            max_iter=1,
        )


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
        SemiexperimentalFitRequest(xyz, (observation,), fixed_parameters=("R(1,2) Frozen",)),
        max_iter=1,
    )

    active_count = sum(parameter.active for parameter in result.parameters)
    assert any("R(1,2)" in parameter.name and parameter.active for parameter in result.parameters)
    assert result.jacobian.shape[1] < active_count


@pytest.mark.parametrize("coordinate_model", ("gic", "cartesian_symmetry"))
def test_semiexperimental_primitive_constraints_preserve_fixed_values(tmp_path, coordinate_model):
    xyz = tmp_path / f"water_{coordinate_model}.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.757000 0.586000",
                "H 0.000000 -0.757000 0.586000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    target_coords = coords.copy()
    target_coords[1:] *= 1.08
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, target_coords))),
    )

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(
            xyz,
            (observation,),
            fixed_parameters=("R(1,2) Frozen",),
            coordinate_model=coordinate_model,
        ),
        max_iter=8,
        outdir=tmp_path / f"fixed_{coordinate_model}",
    )

    initial_r12 = float(np.linalg.norm(coords[0] - coords[1]))
    final = result.final_coordinates_angstrom
    assert float(np.linalg.norm(final[0] - final[1])) == pytest.approx(initial_r12, abs=1.0e-8)
    assert float(np.linalg.norm(final[0] - final[2])) == pytest.approx(initial_r12, abs=1.0e-8)


@pytest.mark.parametrize("coordinate_model", ("gic", "cartesian_symmetry"))
def test_semiexperimental_linear_primitive_constraints_preserve_combinations(tmp_path, coordinate_model):
    xyz = tmp_path / f"water_linear_{coordinate_model}.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.757000 0.586000",
                "H 0.000000 -0.757000 0.586000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    target_coords = coords.copy()
    target_coords[1] *= 1.12
    target_coords[2] *= 0.95
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, target_coords))),
    )

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(
            xyz,
            (observation,),
            fixed_parameters=("DR(Frozen,Value=0.0)=R[1,3]-R[1,2]",),
            coordinate_model=coordinate_model,
        ),
        max_iter=10,
        outdir=tmp_path / f"linear_fixed_{coordinate_model}",
    )

    final = result.final_coordinates_angstrom
    r12 = float(np.linalg.norm(final[0] - final[1]))
    r13 = float(np.linalg.norm(final[0] - final[2]))
    assert r13 - r12 == pytest.approx(0.0, abs=1.0e-8)


@pytest.mark.parametrize("coordinate_model", ("gic", "cartesian_symmetry"))
def test_semiexperimental_gaussian_expression_constraints_preserve_combinations(tmp_path, coordinate_model):
    xyz = tmp_path / f"water_gic_expression_{coordinate_model}.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.757000 0.586000",
                "H 0.000000 -0.757000 0.586000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    atoms = ["O", "H", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    target_coords = coords.copy()
    target_coords[1] *= 1.12
    target_coords[2] *= 0.95
    observation = IsotopologueObservation(
        "parent",
        RotationalConstants(*rotational_constants_MHz(_structure(atoms, target_coords))),
    )

    result = fit_semiexperimental_geometry(
        SemiexperimentalFitRequest(
            xyz,
            (observation,),
            fixed_parameters=("DR(Frozen,Value=0.0)=R[1,3]-R[1,2]",),
            coordinate_model=coordinate_model,
        ),
        max_iter=10,
        outdir=tmp_path / f"gic_expression_fixed_{coordinate_model}",
    )

    final = result.final_coordinates_angstrom
    r12 = float(np.linalg.norm(final[0] - final[1]))
    r13 = float(np.linalg.norm(final[0] - final[2]))
    assert r13 - r12 == pytest.approx(0.0, abs=1.0e-8)


def test_semiexperimental_gaussian_expression_constraints_can_reference_gic_names(tmp_path):
    atoms = ("O", "H", "H")
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    prims, u_matrix, labels = _gic_model(
        coords,
        z_numbers,
        backend=_make_gicforge_backend(atoms, tmp_path),
    )

    constraints = _gic_expression_constraints_from_patterns(("QFIX=[GIC001+2*GIC002] Value=0.0",))
    values = _gic_expression_constraint_values(constraints, coords, prims, u_matrix, labels)
    targets = _gic_expression_constraint_targets(constraints, coords, prims, u_matrix, labels)

    assert len(constraints) == 1
    assert constraints[0].name == "QFIX"
    assert targets[0] == pytest.approx(0.0)
    assert values[0] == pytest.approx(_gic_expression_constraint_values(constraints, coords, prims, u_matrix, labels)[0])


def test_semiexperimental_gaussian_gic_keyword_syntax_is_accepted():
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    prims = ()
    u_matrix = np.zeros((0, 0))
    labels = ()
    constraints = _gic_expression_constraints_from_patterns(
        (
            "HOH(Frozen)=A(2,1,3)",
            "QFIX(Frozen,Value=0.0)=R[1,3]-R[1,2]",
            "R[1,3]-R[1,2] Freeze",
            "BondFix(F,Value=1.09)=B[1,2]",
            "AngleFix(Frozen,Value=90.0)=Angle(2,1,3)",
            "CartZ(Frozen,Value=1.0)=Z(2)",
            "DD(Frozen,Value=0.0)=DotDiff(1,2,3,4)",
        )
    )
    values = _gic_expression_constraint_values(constraints, coords, prims, u_matrix, labels)
    targets = _gic_expression_constraint_targets(constraints, coords, prims, u_matrix, labels)

    assert _fixed_primitives_from_patterns(("HOH(Frozen)=A(2,1,3)",))[0].kind == "angle"
    assert [item.name for item in constraints] == ["QFIX", "R[1,3]-R[1,2] Freeze", "BondFix", "AngleFix", "CartZ", "DD"]
    assert targets[0] == pytest.approx(0.0)
    assert targets[1] == pytest.approx(values[1])
    assert targets[2] == pytest.approx(1.09)
    assert targets[3] == pytest.approx(np.pi / 2.0)
    assert targets[4] == pytest.approx(1.0)
    assert targets[5] == pytest.approx(0.0)
    inactive = _gic_expression_constraints_from_patterns(("RPck001(Inactive,Value=1.0)=D(1,2,3,4)",))
    assert inactive == ()
    assert _gic_fixed_patterns(("RPck001(Inactive,Value=1.0)=D(1,2,3,4)",)) == ()


def test_semiexperimental_constraint_jacobians_match_finite_differences():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.1, 0.0],
            [1.0, 1.0, 0.2],
        ],
        dtype=float,
    )
    fixed = _fixed_primitives_from_patterns(("R(1,2) Frozen",))
    fixed_targets = np.array([1.0], dtype=float)
    linear = _linear_primitive_constraints_from_patterns(("Linear(bond(1,3)-bond(1,2)=0.0)",))
    expressions = _gic_expression_constraints_from_patterns(
        (
            "DR(Frozen,Value=0.0)=R(1,3)-R(1,2)",
            "ZLOCK(Frozen,Value=0.2)=Z(4)",
        )
    )
    targets = _gic_expression_constraint_targets(expressions, coords, (), np.zeros((0, 0)), ())

    analytic = _combined_primitive_constraint_b_matrix(
        coords,
        fixed,
        linear,
        expression_constraints=expressions,
        prims=(),
        u_matrix=np.zeros((0, 0)),
        labels=(),
    )
    numeric = _finite_difference_constraint_b_matrix(
        coords,
        fixed,
        fixed_targets,
        linear,
        expression_constraints=expressions,
        expression_targets=targets,
        prims=(),
        u_matrix=np.zeros((0, 0)),
        labels=(),
    )

    assert analytic.shape == numeric.shape
    assert analytic == pytest.approx(numeric, abs=2.0e-6)


def test_semiexperimental_public_constraint_and_mass_apis():
    atoms = ("O", "H", "H")
    observation = IsotopologueObservation(
        "D1",
        RotationalConstants(1.0, 2.0, 3.0),
        substitutions={2: 2},
    )
    fixed, linear, expressions, definitions = parse_gaussian_style_constraints(
        (
            "R(1,2) Frozen",
            "Linear(bond(1,3)-bond(1,2)=0.0)",
            "R12=R(1,2)",
            "LOCK(Frozen,Value=1.0)=R12",
        )
    )

    assert len(fixed) == 1
    assert len(linear) == 1
    assert len(expressions) == 1
    assert len(definitions) == 2
    assert mass_vector_for_observation(atoms, observation)[1] > mass_vector_for_observation(atoms, observation)[2]
    assert callable(finite_difference_constraint_b_matrix)


def test_semiexperimental_uncertainty_diagnostics_report_cutoff_sensitivity():
    labels = ("q1", "q2")
    jac = np.array([[1.0, 0.0], [0.0, 1.0e-9], [1.0, 1.0e-9]], dtype=float)
    residual = np.array([0.1, -0.2, 0.1], dtype=float)

    text = _uncertainty_diagnostics_csv(labels, jac, residual)

    assert "relative_cutoff" in text
    assert "sigma_ratio_to_default" in text
    assert "rel_1e-8" in text


def test_semiexperimental_gaussian_gic_definitions_are_reusable_constraints():
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    prims = ()
    u_matrix = np.zeros((0, 0))
    labels = ()
    records = (
        "R12=R(1,2)",
        "R13=R(1,3)",
        "ADEF=A(2,1,3)",
        "DR(Frozen,Value=0.0)=R13-R12",
        "COSANG(Frozen,Value=0.0)=cos(A(2,1,3))",
        "AFIX(Frozen,Value=90.0)=ADEF",
    )

    definitions = _gic_expression_definitions_from_patterns(records)
    constraints = _gic_expression_constraints_from_patterns(records)
    values = _gic_expression_constraint_values(
        constraints,
        coords,
        prims,
        u_matrix,
        labels,
        definitions=definitions,
    )
    targets = _gic_expression_constraint_targets(
        constraints,
        coords,
        prims,
        u_matrix,
        labels,
        definitions=definitions,
    )

    assert [definition.name for definition in definitions] == ["R12", "R13", "ADEF", "DR", "COSANG", "AFIX"]
    assert [constraint.name for constraint in constraints] == ["DR", "COSANG", "AFIX"]
    assert values[0] == pytest.approx(0.0)
    assert values[1] == pytest.approx(0.0, abs=1.0e-12)
    assert targets == pytest.approx([0.0, 0.0, np.pi / 2.0])


def test_semiexperimental_gaussian_gic_definition_cycles_are_rejected():
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    records = ("ADEF=BDEF+1", "BDEF=ADEF+1", "LOCK(Frozen,Value=0.0)=ADEF")
    definitions = _gic_expression_definitions_from_patterns(records)
    constraints = _gic_expression_constraints_from_patterns(records)

    with pytest.raises(ValueError, match="Cyclic GIC expression definition"):
        _gic_expression_constraint_values(
            constraints,
            coords,
            (),
            np.zeros((0, 0)),
            (),
            definitions=definitions,
        )


def test_semiexperimental_value_constraints_do_not_use_starting_geometry(tmp_path):
    atoms = ["O", "H", "H"]
    base = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    starts = [
        np.vstack((base[0], 0.94 * base[1], 0.94 * base[2])),
        np.vstack((base[0], 1.06 * base[1], 1.06 * base[2])),
    ]
    target = 0.9572
    finals = []
    for idx, coords in enumerate(starts):
        xyz = tmp_path / f"water_value_{idx}.xyz"
        xyz.write_text(
            "\n".join(
                [
                    "3",
                    "water",
                    *(
                        f"{atom} {row[0]:.8f} {row[1]:.8f} {row[2]:.8f}"
                        for atom, row in zip(atoms, coords)
                    ),
                ]
            )
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
                fixed_parameters=(
                    f"R12=[R(1,2)] Value={target}",
                    f"R13=[R(1,3)] Value={target}",
                ),
            ),
            max_iter=0,
        )
        finals.append(result.final_coordinates_angstrom)

    for final in finals:
        assert float(np.linalg.norm(final[0] - final[1])) == pytest.approx(target, abs=1.0e-8)
        assert float(np.linalg.norm(final[0] - final[2])) == pytest.approx(target, abs=1.0e-8)


def test_semiexperimental_primitive_constraints_expand_by_symmetry(tmp_path):
    atoms = ("O", "H", "H")
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.7570, 0.5860], [0.0, -0.7570, 0.5860]])
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    prims, _u_matrix, _labels = _gic_model(
        coords,
        z_numbers,
        backend=_make_gicforge_backend(atoms, tmp_path),
    )
    fixed = _fixed_primitives_from_patterns(("R(1,2) Frozen",))

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
                ParameterClassConstraint("OH_stretches", ("R(1,2)", "R(1,3)"), "shared"),
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
    assert "Warnings" in report_text
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
        parameter_classes=(ParameterClassConstraint("OH_stretches", ("R(1,2)", "R(1,3)"), "shared"),),
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
    assert any("D(" in label for label in preview.gic_labels)


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
