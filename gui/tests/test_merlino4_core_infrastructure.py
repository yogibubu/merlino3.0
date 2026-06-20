from __future__ import annotations

import json

from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_core import build_run_manifest, ensure_workspace, load_config, write_default_config
from merlino_core.cli import main as merlino_cli
from merlino_gaussian import summarize_gaussian_log
from merlino_semiexp import IsotopologueObservation, RotationalConstants, write_observations_csv
from merlino_gui import discover_manifests


def test_workspace_config_and_manifest_contracts(tmp_path):
    layout = ensure_workspace(tmp_path / "project")
    config_path = write_default_config(layout.root / "merlino.toml")
    input_file = layout.inputs / "input.txt"
    output_file = layout.outputs / "output.txt"
    input_file.write_text("input\n", encoding="utf-8")
    output_file.write_text("output\n", encoding="utf-8")

    config = load_config(config_path)
    manifest_path = build_run_manifest(
        workflow="demo",
        status="completed",
        run_dir=layout.new_run_dir("demo"),
        inputs={"input": input_file},
        outputs={"output": output_file},
        parameters={"answer": 42},
        backend={"name": "python"},
    ).write()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert config.gaussian_executable == "gdv"
    assert layout.inputs.exists()
    assert layout.reports.exists()
    assert manifest["schema_version"] == "merlino.run.v1"
    assert manifest["workflow"] == "demo"
    assert "input" in manifest["input_sha256"]
    assert "output" in manifest["output_sha256"]


def test_gaussian_log_summary_parser(tmp_path):
    log = tmp_path / "test.log"
    log.write_text(
        "\n".join(
            [
                " Standard orientation:",
                " ---------------------------------------------------------------------",
                " Center     Atomic      Atomic             Coordinates (Angstroms)",
                " Number     Number       Type             X           Y           Z",
                " ---------------------------------------------------------------------",
                "      1          8           0        0.000000    0.000000    0.000000",
                "      2          1           0        0.000000    0.000000    0.960000",
                " ---------------------------------------------------------------------",
                " SCF Done:  E(RB3LYP) = -76.123456 A.U.",
                " Frequencies -- 1000.0 1500.0 2000.0",
                " QPck001 PhiP001 RPck001",
                " Normal termination of Gaussian 16",
            ]
        ),
        encoding="utf-8",
    )

    summary = summarize_gaussian_log(log)

    assert summary.normal_termination
    assert summary.scf_energies_hartree == (-76.123456,)
    assert summary.standard_orientation_count == 1
    assert summary.puckering_marker_count == 3
    assert summary.frequencies_cm == (1000.0, 1500.0, 2000.0)
    assert len(summary.last_orientation) == 2


def test_merlino_cli_init_vci_and_dvr_args(tmp_path):
    project = tmp_path / "project"
    assert merlino_cli(["init", str(project)]) == 0
    assert (project / "merlino.toml").exists()
    assert (project / "runs").exists()

    qff = tmp_path / "field.qff"
    qff.write_text(
        "\n".join(["FREQ 1 1000.0", "FREQ 2 1500.0", "QUARTIC 1 1 1 1 0.8"]) + "\n",
        encoding="utf-8",
    )
    report = tmp_path / "vci_report.txt"
    csv_dir = tmp_path / "vci_csv"
    run_dir = tmp_path / "vci_run"
    assert merlino_cli([
        "vci",
        "--qff",
        str(qff),
        "--max-quanta",
        "2",
        "--roots",
        "3",
        "--out",
        str(report),
        "--run-dir",
        str(run_dir),
        "--csv-dir",
        str(csv_dir),
    ]) == 0
    assert "VPT2/VCI comparison" in report.read_text(encoding="utf-8")
    assert (run_dir / "vpt2_vci_manifest.json").exists()
    assert (csv_dir / "vpt2_vci_comparison.csv").exists()

    log = tmp_path / "scan.log"
    log.write_text("Normal termination of Gaussian 16\n", encoding="utf-8")
    outdir = tmp_path / "dvr_out"
    figdir = tmp_path / "dvr_fig"
    assert merlino_cli([
        "dvr-args",
        "--repo-root",
        str(tmp_path),
        "--log",
        str(log),
        "--outdir",
        str(outdir),
        "--figdir",
        str(figdir),
        "--prefix",
        "scan",
    ]) == 0
    assert (outdir / "scan_manifest.json").exists()


def test_merlino_cli_semiexp(tmp_path):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.987200",
                "H 0.906600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target = Structure.from_atoms_coords(
        ["O", "H", "H"],
        [(0.0, 0.0, 0.0), (0.0, 0.0, 0.9572), (0.9266, 0.0, -0.2396)],
    )
    observations = (
        IsotopologueObservation("parent", RotationalConstants(*rotational_constants_MHz(target))),
    )
    obs_csv = write_observations_csv(tmp_path / "observations.csv", observations)
    outdir = tmp_path / "semiexp"

    assert merlino_cli(["semiexp", "--xyz", str(xyz), "--observations", str(obs_csv), "--outdir", str(outdir)]) == 0

    assert (outdir / "semiexp_geometry.xyz").exists()
    assert (outdir / "semiexp_parameters.csv").exists()
    assert (outdir / "semiexp_residuals.csv").exists()
    assert (outdir / "semiexp_manifest.json").exists()


def test_merlino_cli_gic_gaussian_summary_and_backends(tmp_path):
    executable = tmp_path / "fake_gicforge.sh"
    executable.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -e",
                "printf 'legacy report\\n' > provout",
                "printf 'gaussian input\\n' > gauin",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    workdir = tmp_path / "gic"
    workdir.mkdir()
    (workdir / "provin").write_text("input\n", encoding="utf-8")

    assert merlino_cli(["gic", "--workdir", str(workdir), "--executable", str(executable)]) == 0
    assert (workdir / "gicforge_manifest.json").exists()
    assert (workdir / "gicforge.out").exists()

    log = tmp_path / "gaussian.log"
    log.write_text(
        "SCF Done:  E(RHF) = -1.000000 A.U.\nNormal termination of Gaussian 16\n",
        encoding="utf-8",
    )
    assert merlino_cli(["gaussian-summary", str(log)]) == 0
    assert merlino_cli(["backends"]) == 0
    assert merlino_cli(["compare-backends"]) == 0


def test_manifest_discovery(tmp_path):
    run_dir = tmp_path / "runs" / "demo"
    run_dir.mkdir(parents=True)
    manifest = build_run_manifest(workflow="demo", status="completed", run_dir=run_dir)
    path = manifest.write()

    entries = discover_manifests(tmp_path)

    assert entries[0].path == path
    assert entries[0].workflow == "demo"
