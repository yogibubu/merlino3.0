from __future__ import annotations

from pathlib import Path

import numpy as np

from merlino_semiexp import (
    EnsembleClassCorrection,
    EnsembleMolecule,
    IsotopologueObservation,
    RotationalConstants,
    SemiexperimentalFitRequest,
    fit_ensemble_job,
    fit_ensemble_class_corrections,
    run_ensemble_leave_one_molecule_out,
    run_ensemble_prior_comparison,
    run_ensemble_prior_scan,
    run_ensemble_synthon_threshold_scan,
    suggest_ensemble_class_priors,
)
from merlino_semiexp.fit import _rotational_constants_for_substitution
from merlino_semiexp.ensemble import _class_projection, _primitive_signature


def test_ensemble_class_correction_reduces_shared_residual(tmp_path):
    water = tmp_path / "water.xyz"
    water.write_text(
        "\n".join(
            [
                "3",
                "water-like reference",
                "O 0.000000 0.000000 0.000000",
                "H 0.958400 0.000000 0.000000",
                "H -0.239000 0.927000 0.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sulfide = tmp_path / "sulfide.xyz"
    sulfide.write_text(
        "\n".join(
            [
                "3",
                "sulfide-like reference",
                "S 0.000000 0.000000 0.000000",
                "H 1.330000 0.000000 0.000000",
                "H -0.335000 1.180000 0.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    water_constants = _rotational_constants_for_substitution(
        ("O", "H", "H"),
        np.array([[0.0, 0.0, 0.0], [0.9584, 0.0, 0.0], [-0.239, 0.927, 0.0]], dtype=float),
        {},
    )
    sulfide_constants = _rotational_constants_for_substitution(
        ("S", "H", "H"),
        np.array([[0.0, 0.0, 0.0], [1.33, 0.0, 0.0], [-0.335, 1.18, 0.0]], dtype=float),
        {},
    )

    # A common stretch-class bias lowers the large rotational components for
    # both references.  The exact value is unimportant; the ensemble layer
    # should find one shared correction that improves both residual blocks.
    water_obs = water_constants + np.array([-1200.0, -400.0, 0.0])
    sulfide_obs = sulfide_constants + np.array([-180.0, -70.0, 0.0])

    molecules = (
        EnsembleMolecule(
            "water",
            SemiexperimentalFitRequest(
                water,
                (_observation("parent", water_obs),),
                observable="rotational_constants",
                rotational_components="AB",
            ),
        ),
        EnsembleMolecule(
            "sulfide",
            SemiexperimentalFitRequest(
                sulfide,
                (_observation("parent", sulfide_obs),),
                observable="rotational_constants",
                rotational_components="AB",
            ),
        ),
    )
    result = fit_ensemble_class_corrections(
        molecules,
        (EnsembleClassCorrection("XH_stretch_bias", ("R(",)),),
        outdir=tmp_path / "ensemble",
    )

    assert result.rank == 1
    assert result.condition_number < np.inf
    assert result.acceptance.status == "accepted"
    assert result.acceptance.accepted is True
    assert result.diagnostics.rank == result.rank
    assert result.diagnostics.n_columns == 1
    assert result.diagnostics.residual_degrees_of_freedom >= 0
    assert result.weighted_rms_after < result.weighted_rms_before
    assert result.molecule_blocks[0].matched_counts["XH_stretch_bias"] == 1
    assert (tmp_path / "ensemble" / "ensemble_class_corrections.csv").exists()
    manifest = (tmp_path / "ensemble" / "ensemble_manifest.json")
    assert manifest.exists()
    manifest_text = manifest.read_text(encoding="utf-8")
    assert "numerical_diagnostics" in manifest_text
    assert "q_SE(m,k) = q_QC(m,k) + Delta[class(m,k)]" in manifest_text
    assert "high_correlation_reject_threshold" in manifest_text

    prior_result = fit_ensemble_class_corrections(
        molecules,
        (EnsembleClassCorrection("XH_stretch_bias", ("R(",), prior_value=0.0, prior_sigma=1.0e-4),),
    )
    assert abs(prior_result.corrections["XH_stretch_bias"]) < abs(result.corrections["XH_stretch_bias"])
    assert "XH_stretch_bias" in prior_result.prior_residual_after


def test_anhydrides_ensemble_job_is_rank_consistent(tmp_path):
    root = Path(__file__).resolve().parents[2]
    job = root / "examples" / "semiexp" / "anhydrides_ensemble" / "anhydrides_ensemble.mse-ensemble.toml"
    result = fit_ensemble_job(job, outdir=tmp_path / "anhydrides")

    assert result.rank == 8
    assert result.condition_number < 1.0e4
    assert result.acceptance.status == "review"
    assert result.acceptance.accepted is True
    assert any("CO_carbonyl/CO_single" in item for item in result.acceptance.review_items)
    assert result.weighted_rms_after < result.weighted_rms_before
    assert {item.name for item in result.classes} == {
        "CC_short",
        "CC_long",
        "CO_carbonyl",
        "CO_single",
        "CCC_bend",
        "CCO_bend",
        "COC_bend",
        "OCO_bend",
    }
    assert {"CCC_bend", "CCO_bend", "COC_bend", "OCO_bend"} <= set(result.prior_residual_after)
    assert {block.molecule for block in result.molecule_blocks} == {
        "maleic_anhydride",
        "phthalic_anhydride",
        "succinic_anhydride",
    }
    assert (tmp_path / "anhydrides" / "ensemble_class_corrections.txt").exists()
    assert (tmp_path / "anhydrides" / "ensemble_class_report.csv").exists()


def test_anhydrides_parent_only_ensemble_is_rank_consistent(tmp_path):
    root = Path(__file__).resolve().parents[2]
    job = root / "examples" / "semiexp" / "anhydrides_parent_only" / "anhydrides_parent_only.mse-ensemble.toml"
    result = fit_ensemble_job(job, outdir=tmp_path / "anhydrides_parent_only")

    assert result.rank == 4
    assert len(result.classes) == 4
    assert result.condition_number < 1.0e3
    assert result.acceptance.status == "review"
    assert result.acceptance.accepted is True
    assert any("CO_carbonyl/CO_single" in item for item in result.acceptance.review_items)
    assert result.weighted_rms_after < result.weighted_rms_before
    assert {block.molecule for block in result.molecule_blocks} == {
        "maleic_anhydride",
        "phthalic_anhydride",
        "succinic_anhydride",
    }
    assert (tmp_path / "anhydrides_parent_only" / "ensemble_class_report.csv").exists()


def test_anhydrides_prior_comparison_writes_variants(tmp_path):
    root = Path(__file__).resolve().parents[2]
    job = root / "examples" / "semiexp" / "anhydrides_ensemble" / "anhydrides_ensemble.mse-ensemble.toml"
    results = run_ensemble_prior_comparison(job, tmp_path / "comparison")

    assert set(results) == {"no_prior", "soft_prior", "hard_constraint"}
    assert results["no_prior"].condition_number > results["soft_prior"].condition_number
    assert results["no_prior"].acceptance.status == "rejected"
    assert results["soft_prior"].acceptance.status == "review"
    assert results["hard_constraint"].acceptance.status == "review"
    assert results["soft_prior"].rank == 8
    assert results["hard_constraint"].rank == 4
    assert (tmp_path / "comparison" / "ensemble_prior_comparison.csv").exists()
    assert (tmp_path / "comparison" / "soft_prior" / "ensemble_manifest.json").exists()
    assert (tmp_path / "comparison" / "prior_scan" / "prior_sigma_scan.csv").exists()
    assert (tmp_path / "comparison" / "leave_one_molecule_out" / "leave_one_molecule_out.csv").exists()

    scan = run_ensemble_prior_scan(job, tmp_path / "scan", sigmas=(1.0e-4, 1.0e-3, 1.0e-2))
    assert len(scan) == 3
    loo = run_ensemble_leave_one_molecule_out(job, tmp_path / "loo")
    assert {row["heldout"] for row in loo} == {"maleic_anhydride", "phthalic_anhydride", "succinic_anhydride"}
    assert all(float(row["transferability_score"]) > 0.0 for row in loo)


def test_glycine_conformer_ensemble_uses_continuous_synthon_atom_types(tmp_path):
    root = Path(__file__).resolve().parents[2]
    job = root / "examples" / "semiexp" / "glycine_ensemble" / "glycine_conformers_synthon.mse-ensemble.toml"
    result = fit_ensemble_job(job, outdir=tmp_path / "glycine_synthon")

    assert result.rank == 10
    assert result.condition_number < 1.0e4
    assert result.acceptance.status == "accepted"
    assert result.acceptance.accepted is True
    assert result.weighted_rms_after < result.weighted_rms_before
    assert {block.molecule for block in result.molecule_blocks} == {"glycine_I", "glycine_II"}
    assert all(block.matched_counts["Ccarb_Ocarbonyl"] > 0 for block in result.molecule_blocks)
    assert all(block.matched_counts["Ccarb_Ohydroxyl"] > 0 for block in result.molecule_blocks)
    assert (tmp_path / "glycine_synthon" / "ensemble_class_report.csv").exists()

    scan = run_ensemble_synthon_threshold_scan(job, tmp_path / "glycine_scan", thresholds=(0.01, 0.035, 0.1))
    assert [row["status"] for row in scan] == ["ok", "ok", "ok"]
    assert [row["acceptance_status"] for row in scan] == ["accepted", "accepted", "rejected"]
    assert int(scan[0]["rank"]) == int(scan[1]["rank"]) == 10
    assert int(scan[2]["rank"]) < 10
    assert (tmp_path / "glycine_scan" / "synthon_threshold_scan.csv").exists()


def test_auto_prior_policy_only_regularizes_angular_classes():
    classes = (
        EnsembleClassCorrection("CC", kind="stretch", atom_symbols=("C", "C")),
        EnsembleClassCorrection("CCC", kind="bend", atom_symbols=("C", "C", "C")),
        EnsembleClassCorrection("CCT", kind="torsion", atom_symbols=("C", "C")),
    )
    suggested = suggest_ensemble_class_priors(classes)
    assert suggested[0].prior_sigma is None
    assert suggested[1].prior_sigma == 1.0e-3
    assert suggested[2].prior_sigma == 2.0e-3


def test_ensemble_coordinate_signatures_preserve_chemical_centers():
    atoms = ("C", "N", "C", "H", "O")

    assert _primitive_signature("torsion", (4, 1, 2, 3), atoms) == ("C", "N")
    assert _primitive_signature("torsion", (3, 2, 1, 4), atoms) == ("C", "N")
    assert _class_projection(
        EnsembleClassCorrection("CN_torsion", kind="torsion", atom_symbols=("N", "C")),
        "GIC [ 0.5*D(  4,  1,  2,  3)-0.5*D(  3,  2,  1,  4)]",
        atoms,
    ) == 0.0

    assert _primitive_signature("out_of_plane", (2, 1, 3, 4), atoms) == ("N",)
    assert _class_projection(
        EnsembleClassCorrection("N_oop", kind="out_of_plane", atom_symbols=("N",)),
        "GIC [ 0.75*U(  2,  1,  3,  4)+0.25*U(  1,  2,  3,  5)]",
        atoms,
    ) == 0.75


def test_ensemble_value_and_synthon_selectors_refine_classes():
    atoms = ("C", "O", "C")
    coords = np.array([[0.0, 0.0, 0.0], [1.21, 0.0, 0.0], [2.63, 0.0, 0.0]], dtype=float)
    label = "GIC [ 1.0*R(  1,  2)+1.0*R(  2,  3)]"

    assert _class_projection(
        EnsembleClassCorrection("CO_carbonyl", kind="stretch", atom_symbols=("C", "O"), value_max=1.28),
        label,
        atoms,
        coords=coords,
    ) == 1.0
    assert _class_projection(
        EnsembleClassCorrection("CO_single", kind="stretch", atom_symbols=("C", "O"), value_min=1.28),
        label,
        atoms,
        coords=coords,
    ) == 1.0
    assert _class_projection(
        EnsembleClassCorrection("synthon_filtered", kind="stretch", atom_symbols=("C", "O"), synthon_signatures=("6-0-3-3",)),
        label,
        atoms,
        coords=coords,
        synthon_signatures=("6-0-3-3", "8-0-3-2", "6-0-4-4"),
    ) == 1.0
    assert _class_projection(
        EnsembleClassCorrection(
            "zeff_filtered",
            kind="stretch",
            atom_symbols=("C", "O"),
            synthon_zeff=(6.35, 8.20),
            synthon_threshold=0.10,
        ),
        label,
        atoms,
        coords=coords,
        synthon_zeff=(6.33, 8.23, 6.48),
    ) == 1.0
    assert _class_projection(
        EnsembleClassCorrection(
            "zeff_rejected",
            kind="stretch",
            atom_symbols=("C", "O"),
            synthon_zeff=(6.80, 8.20),
            synthon_threshold=0.10,
        ),
        label,
        atoms,
        coords=coords,
        synthon_zeff=(6.33, 8.23, 6.48),
    ) == 0.0


def test_ensemble_job_parser_normalizes_synthon_signature_strings(tmp_path):
    job = tmp_path / "ensemble.mse-ensemble.toml"
    geom = tmp_path / "water.xyz"
    geom.write_text(
        "3\nwater\nO 0 0 0\nH 0.95 0 0\nH -0.2 0.9 0\n",
        encoding="utf-8",
    )
    molecule_job = tmp_path / "water.mse.toml"
    molecule_job.write_text(
        f"""
schema = "merlino.semiexp.job.v1"
title = "water"

[files]
geometry = "{geom}"

[fit]
observable = "rotational_constants"
rotational_components = "AB"

[[isotopologues]]
label = "parent"

[isotopologues.constants]
A_MHz = 1.0
B_MHz = 1.0
C_MHz = 1.0
""",
        encoding="utf-8",
    )
    job.write_text(
        f"""
schema = "merlino.semiexp.ensemble.v1"
title = "parser"

[acceptance]
require_full_rank = false
max_condition_number = 1000000.0
min_residual_degrees_of_freedom = 2
min_molecule_support = 1
high_correlation_review_threshold = 0.95
high_correlation_reject_threshold = 0.999

[[molecules]]
name = "water"
job = "{molecule_job}"

[[classes]]
name = "synthon"
kind = "stretch"
atoms = "O,H"
synthon_signatures = "8-0-2-0|1-0-1-0"
""",
        encoding="utf-8",
    )

    from merlino_semiexp.ensemble import read_ensemble_job

    parsed = read_ensemble_job(job)
    assert parsed.classes[0].synthon_signatures == ("8-0-2-0", "1-0-1-0")
    assert parsed.acceptance_policy.require_full_rank is False
    assert parsed.acceptance_policy.max_condition_number == 1.0e6
    assert parsed.acceptance_policy.min_residual_degrees_of_freedom == 2
    assert parsed.acceptance_policy.min_molecule_support == 1
    assert parsed.acceptance_policy.high_correlation_review_threshold == 0.95
    assert parsed.acceptance_policy.high_correlation_reject_threshold == 0.999


def test_ensemble_cli_writes_canonical_run_manifest(tmp_path):
    from merlino_core.cli import main

    root = Path(__file__).resolve().parents[2]
    job = root / "examples" / "semiexp" / "anhydrides_ensemble" / "anhydrides_ensemble.mse-ensemble.toml"
    outdir = tmp_path / "cli_anhydrides"
    status = main(["semiexp-ensemble", "--job", str(job), "--outdir", str(outdir)])

    assert status == 0
    run_manifest = outdir / "run_manifest.json"
    assert run_manifest.exists()
    text = run_manifest.read_text(encoding="utf-8")
    assert '"workflow": "semiexp_ensemble"' in text
    assert '"scientific_manifest"' in text
    assert '"input_sha256"' in text
    assert '"output_sha256"' in text


def _observation(label: str, constants: np.ndarray) -> IsotopologueObservation:
    return IsotopologueObservation(
        label,
        RotationalConstants(float(constants[0]), float(constants[1]), float(constants[2])),
    )
