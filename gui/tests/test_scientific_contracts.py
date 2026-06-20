from __future__ import annotations

import pytest

from merlino_semiexp import (
    CorrectedRotationalConstants,
    IsotopologueObservation,
    RotationalConstants,
    SemiexperimentalFitRequest,
    VibrationalCorrection,
    corrected_constants_rows,
    parse_substitutions,
    read_observations_csv,
    write_observations_csv,
)
from merlino_core import repo_root
from merlino_vpt2_vci import DavidsonSettings, ForceFieldSource, VCIRequest, inventory_vpt2_vci_backends


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
    assert inventory.harmonic_internal_source is not None
    assert inventory.harmonic_internal_source.name == "gf.f"
    assert inventory.active_fortran_sources == ()
    assert inventory.davidson_backend is None
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
    assert rows[0]["A_e_MHz"] == 999.0
    assert rows[1]["C_e_MHz"] == 587.5


def test_semiexperimental_substitution_parser():
    assert parse_substitutions("2:13;5:18") == {2: 13, 5: 18}
    with pytest.raises(ValueError):
        parse_substitutions("0:13")
