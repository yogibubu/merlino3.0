from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from advanced.vpt2_vci_window import VPT2VCIWindow


@pytest.mark.usefixtures("qtbot")
def test_vpt2_vci_window_defaults(tmp_path, qtbot):
    window = VPT2VCIWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Merlino GF / VPT2-VCI"
    assert "gauin.fchk" in window.fchk_edit.text()
    assert window.max_quanta_edit.text() == "2"
    assert window.roots_edit.text() == "6"
    assert window.qff_edit.text() == ""


@pytest.mark.usefixtures("qtbot")
def test_vpt2_vci_window_selects_latest_fchk(tmp_path, qtbot):
    old_fchk = tmp_path / "old.fchk"
    new_fchk = tmp_path / "gauin.fchk"
    old_fchk.write_text("old\n", encoding="utf-8")
    new_fchk.write_text("new\n", encoding="utf-8")
    os.utime(old_fchk, (1, 1))
    os.utime(new_fchk, (2, 2))

    window = VPT2VCIWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window._select_latest_fchk(show_message=False)

    assert window.fchk_edit.text() == str(new_fchk)


@pytest.mark.usefixtures("qtbot")
def test_vpt2_vci_window_runs_gf_on_gaussian_fchk(tmp_path, qtbot):
    fchk = Path("gui/tests/gaussian/h2o.fchk").resolve()
    window = VPT2VCIWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.fchk_edit.setText(str(fchk))

    window.run_gf(show_message=False)

    text = window.output_text.toPlainText()
    assert "GF/PED from Merlino non-redundant GICs" in text
    assert "2169.878" in text
    assert "GIC001" in text
    assert "PED (%)" in text
    manifest = json.loads((tmp_path / "gf_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "merlino.run.v1"
    assert manifest["workflow"] == "gf"

    csv_dir = tmp_path / "gf_csv"
    written = window.export_csvs(csv_dir, show_message=False)
    assert (csv_dir / "gf_frequencies.csv").exists()
    assert "ped.csv" in written


@pytest.mark.usefixtures("qtbot")
def test_vpt2_vci_window_runs_comparison_from_indexed_qff(tmp_path, qtbot):
    qff = tmp_path / "field.qff"
    qff.write_text(
        "\n".join(
            [
                "FREQ 1 1000.0",
                "FREQ 2 1500.0",
                "CUBIC 1 1 2 -2.0",
                "QUARTIC 1 1 1 1 0.8",
                "QUARTIC 2 2 2 2 0.5",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    window = VPT2VCIWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.fchk_edit.setText("")
    window.qff_edit.setText(str(qff))
    window.max_quanta_edit.setText("2")
    window.roots_edit.setText("4")
    window.active_modes_edit.setText("1,2")

    window.run_vpt2_vci(show_message=False)

    text = window.output_text.toPlainText()
    assert "VPT2/VCI comparison" in text
    assert "Modes used in input force field: 2" in text
    assert "VCI basis size" in text
    assert "1000." in text
    manifest = json.loads((tmp_path / "vpt2_vci_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "merlino.run.v1"
    assert manifest["workflow"] == "vpt2_vci"

    csv_dir = tmp_path / "vci_csv"
    written = window.export_csvs(csv_dir, show_message=False)
    assert (csv_dir / "vpt2_vci_comparison.csv").exists()
    assert "comparison.csv" in written


@pytest.mark.usefixtures("qtbot")
def test_vpt2_vci_window_previews_cli(tmp_path, qtbot):
    window = VPT2VCIWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.qff_edit.setText(str(tmp_path / "field.qff"))
    window.max_quanta_edit.setText("3")
    window.roots_edit.setText("5")

    window.preview_cli_command()

    text = window.output_text.toPlainText()
    assert "python -m merlino vci" in text
    assert "--max-quanta 3" in text


@pytest.mark.usefixtures("qtbot")
def test_vpt2_vci_window_exports_report_and_roundtrips_preset(tmp_path, qtbot):
    window = VPT2VCIWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.output_text.setPlainText("sample report")
    report_path = tmp_path / "report.txt"

    written = window.export_report(report_path, show_message=False)

    assert written == report_path
    assert report_path.read_text(encoding="utf-8") == "sample report\n"

    window.max_quanta_edit.setText("4")
    window.roots_edit.setText("5")
    window.active_modes_edit.setText("1,3")
    window.force_threshold_edit.setText("0.5")
    preset_path = tmp_path / "preset.json"

    saved = window.save_preset(preset_path, show_message=False)
    window.max_quanta_edit.setText("1")
    window.roots_edit.setText("1")
    window.active_modes_edit.setText("")
    window.force_threshold_edit.setText("0.0")
    loaded = window.load_preset(preset_path, show_message=False)

    assert saved == preset_path
    assert loaded == preset_path
    assert window.max_quanta_edit.text() == "4"
    assert window.roots_edit.text() == "5"
    assert window.active_modes_edit.text() == "1,3"
    assert window.force_threshold_edit.text() == "0.5"
