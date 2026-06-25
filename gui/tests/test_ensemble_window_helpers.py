from __future__ import annotations

import tomllib
from pathlib import Path

from gui.ensemble_window import _ensemble_job_text, _join_items, _split_items


def test_ensemble_job_writer_keeps_fine_class_selectors_and_absolute_paths(tmp_path):
    source = tmp_path / "jobs" / "ensemble.mse-ensemble.toml"
    source.parent.mkdir()
    data = {
        "title": "ensemble",
        "fit": {"step": 1.0e-4, "rcond": 1.0e-10},
        "acceptance": {
            "require_full_rank": True,
            "max_condition_number": 1.0e7,
            "min_residual_degrees_of_freedom": 2,
            "min_molecule_support": 3,
            "high_correlation_review_threshold": 0.97,
            "high_correlation_reject_threshold": 0.999,
        },
        "molecules": [{"name": "mol", "job": "../mol/mol.mse.toml"}],
    }
    text = _ensemble_job_text(
        data,
        [
            {
                "name": "CO_carbonyl",
                "kind": "stretch",
                "atoms": ["C", "O"],
                "value_max": 1.28,
                "synthon_signatures": ["6-0-3-3"],
                "synthon_zeff": [6.35, 8.20],
                "synthon_threshold": 0.10,
            },
            {
                "name": "CCO_bend",
                "kind": "bend",
                "atoms": ["C", "C", "O"],
                "prior_value": 0.0,
                "prior_sigma": 1.0e-3,
            },
        ],
        source_path=source,
    )

    parsed = tomllib.loads(text)
    assert Path(parsed["molecules"][0]["job"]).is_absolute()
    assert parsed["acceptance"]["require_full_rank"] is True
    assert parsed["acceptance"]["max_condition_number"] == 1.0e7
    assert parsed["acceptance"]["min_residual_degrees_of_freedom"] == 2
    assert parsed["acceptance"]["min_molecule_support"] == 3
    assert parsed["acceptance"]["high_correlation_review_threshold"] == 0.97
    assert parsed["acceptance"]["high_correlation_reject_threshold"] == 0.999
    assert parsed["classes"][0]["value_max"] == 1.28
    assert parsed["classes"][0]["synthon_signatures"] == ["6-0-3-3"]
    assert parsed["classes"][0]["synthon_zeff"] == [6.35, 8.20]
    assert parsed["classes"][0]["synthon_threshold"] == 0.10
    assert parsed["classes"][1]["prior_sigma"] == 1.0e-3


def test_ensemble_editor_item_split_join_roundtrip():
    assert _split_items("C,O;H|N") == ["C", "O", "H", "N"]
    assert _join_items(["C", "O", "H"]) == "C,O,H"
