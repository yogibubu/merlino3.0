from __future__ import annotations

import csv
import ast
import json
import math
from pathlib import Path

import pytest

from merlino_semiexp.paper_benchmarks import SnapshotValidationError, validate_paper_benchmark_snapshot


ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "benchmarks/semiexp_msr/golden/semiexp_paper_regression.json"
GENERATED_SUMMARY = ROOT / "benchmarks/semiexp_msr/generated/paper_benchmark_summary.csv"
GENERATED_PLANAR = ROOT / "benchmarks/semiexp_msr/generated/paper_planar_pair_diagnostics.csv"


def _load_snapshot() -> dict:
    return json.loads(SNAPSHOT.read_text(encoding="utf-8"))


def _run_dir(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _diagnostics(path: Path) -> dict[str, str]:
    with (path / "semiexp_diagnostics.csv").open(newline="", encoding="utf-8") as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle)}


def _rotational_stats(path: Path) -> tuple[int, float, float]:
    with (path / "semiexp_rotational_constants.csv").open(newline="", encoding="utf-8") as handle:
        residuals = [float(row["difference_MHz"]) for row in csv.DictReader(handle)]
    rms = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    return len(residuals), rms, max(abs(value) for value in residuals)


def _minimum_eigenvalue(path: Path) -> float:
    with (path / "semiexp_hessian_eigenvalues.csv").open(newline="", encoding="utf-8") as handle:
        return float(next(csv.DictReader(handle))["eigenvalue"])


def test_paper_regression_snapshot_has_all_benchmarks():
    snapshot = _load_snapshot()
    validate_paper_benchmark_snapshot(snapshot, source=SNAPSHOT)
    assert snapshot["schema"] == "merlino.semiexp.paper_regression.v1"
    assert tuple(snapshot["cases"]) == (
        "glycolaldehyde",
        "glycine_II",
        "cyclopentadiene",
        "nitrobenzene",
        "azulene",
        "norcamphor",
    )
    for case in snapshot["cases"].values():
        assert case["reported_constants"] == case["isotopologues"] * 3
        assert case["rank"] == case["effective_parameters"]
        assert case["minimum_hessian_eigenvalue"] > 0.0
        assert case["rotational_rms_MHz"] > 0.0
        assert case["rotational_max_MHz"] >= case["rotational_rms_MHz"]
        assert not Path(case["run_dir"]).is_absolute()


def test_planar_pair_regression_selects_ab_for_aromatic_planar_cases():
    snapshot = _load_snapshot()
    for system, data in snapshot["planar_pair_diagnostics"].items():
        selected = data["selected"]
        pairs = data["pairs"]
        assert selected == "AB", system
        assert set(pairs) == {"AB", "AC", "BC"}
        best_rms_pair = min(pairs, key=lambda key: pairs[key]["rotational_rms_MHz"])
        assert best_rms_pair == selected
        assert pairs[selected]["rank"] == max(pair["rank"] for pair in pairs.values())
        for pair in pairs.values():
            assert not Path(pair["run_dir"]).is_absolute()


def test_paper_regression_snapshot_validation_rejects_absolute_run_dirs():
    snapshot = _load_snapshot()
    snapshot["cases"]["nitrobenzene"]["run_dir"] = "/tmp/not-portable"

    with pytest.raises(SnapshotValidationError):
        validate_paper_benchmark_snapshot(snapshot)


def test_generated_paper_tables_match_regression_snapshot():
    snapshot = _load_snapshot()
    assert GENERATED_SUMMARY.is_file()
    assert GENERATED_PLANAR.is_file()
    with GENERATED_SUMMARY.open(newline="", encoding="utf-8") as handle:
        rows = {row["system"]: row for row in csv.DictReader(handle)}
    assert tuple(rows) == tuple(snapshot["cases"])
    for name, expected in snapshot["cases"].items():
        row = rows[name]
        assert row["coordinate_model"] == expected["coordinate_model"]
        assert row["rotational_pair"] == expected["rotational_pair"]
        assert int(row["final_gics"]) == expected["final_gics"]
        assert int(row["totally_symmetric_gics"]) == expected["totally_symmetric_gics"]
        assert int(row["effective_parameters"]) == expected["effective_parameters"]
        assert int(row["primitive_constraints"]) == expected["primitive_constraints"]
        assert float(row["rotational_rms_MHz"]) == pytest.approx(expected["rotational_rms_MHz"], rel=1e-12)
        assert float(row["rotational_max_MHz"]) == pytest.approx(expected["rotational_max_MHz"], rel=1e-12)
    with GENERATED_PLANAR.open(newline="", encoding="utf-8") as handle:
        pair_rows = {(row["system"], row["pair"]): row for row in csv.DictReader(handle)}
    for system, diagnostics in snapshot["planar_pair_diagnostics"].items():
        for pair, expected in diagnostics["pairs"].items():
            row = pair_rows[(system, pair)]
            assert row["selected"] == diagnostics["selected"]
            assert int(row["rank"]) == expected["rank"]
            assert float(row["condition_number"]) == pytest.approx(expected["condition_number"], rel=1e-12)
            assert row["run_dir"] == expected["run_dir"]


def test_available_semiexp_outputs_match_paper_regression_snapshot():
    snapshot = _load_snapshot()
    checked = 0
    for name, expected in snapshot["cases"].items():
        run_dir = _run_dir(expected["run_dir"])
        if not run_dir.exists():
            continue
        checked += 1
        diagnostics = _diagnostics(run_dir)
        n_constants, rms, max_abs = _rotational_stats(run_dir)
        assert diagnostics.get("coordinate_model", "gic") == expected["coordinate_model"]
        assert tuple(ast.literal_eval(diagnostics["components"])) == tuple(expected["components"])
        assert int(diagnostics["rank"]) == expected["rank"], name
        assert int(diagnostics["n_optimized_parameters"]) == expected["effective_parameters"], name
        assert int(diagnostics["accepted_steps"]) == expected["accepted_steps"], name
        assert int(diagnostics["rejected_steps"]) == expected["rejected_steps"], name
        assert float(diagnostics["condition_number"]) == pytest.approx(expected["condition_number"], rel=1e-10)
        assert n_constants == expected["reported_constants"], name
        assert rms == pytest.approx(expected["rotational_rms_MHz"], rel=1e-10)
        assert max_abs == pytest.approx(expected["rotational_max_MHz"], rel=1e-10)
        assert _minimum_eigenvalue(run_dir) == pytest.approx(expected["minimum_hessian_eigenvalue"], rel=1e-10)
    assert checked >= 3
