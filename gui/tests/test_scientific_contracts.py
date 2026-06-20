from __future__ import annotations

import pytest

from merlino_core import repo_root
from merlino_semiexp import (
    CorrectedRotationalConstants,
    IsotopologueObservation,
    RotationalConstants,
    SemiexperimentalFitRequest,
    VibrationalCorrection,
)
from merlino_vpt2_vci import DavidsonSettings, ForceFieldSource, VCIRequest, inventory_legacy_fortran


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


def test_vpt2_vci_inventory_records_legacy_vci1d():
    inventory = inventory_legacy_fortran(repo_root(__file__))
    assert inventory.legacy_vci1d is not None
    assert inventory.legacy_vci1d.name == "vci1d.f"
    assert inventory.notes


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
