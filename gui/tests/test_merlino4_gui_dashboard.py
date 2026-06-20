from __future__ import annotations

import pytest
from PySide6.QtCore import Qt

from merlino_gui import DashboardWindow, default_workflows
from merlino_gui.app import build_parser
from merlino_gui.dashboard import workflow_detail_text, workflow_state_lines


def test_default_workflows_include_new_scientific_areas():
    workflows = {item.workflow_id: item for item in default_workflows()}
    assert "vpt2_vci" in workflows
    assert "semiexp_geometry" in workflows
    assert workflows["vpt2_vci"].service == "merlino_vpt2_vci"
    assert workflows["semiexp_geometry"].service == "merlino_semiexp"
    assert "fortran77" in workflows["semiexp_geometry"].backends
    assert workflows["semiexp_geometry"].status == "standard solver"
    assert workflows["gic"].default_backend == "fortran77"


def test_workflow_detail_text_lists_contract_fields():
    workflow = default_workflows()[0]
    text = workflow_detail_text(workflow)
    assert workflow.title in text
    assert "Category:" in text
    assert "Inputs:" in text
    assert "Outputs:" in text
    assert workflow.service in text
    assert "Backend:" in text


def test_workflow_state_reports_manifest_outputs(tmp_path):
    run = tmp_path / "semiexp"
    run.mkdir()
    report = run / "semiexp_report.html"
    report.write_text("<html></html>\n", encoding="utf-8")
    manifest = run / "semiexp_manifest.json"
    manifest.write_text(
        """
{
  "workflow": "semiexperimental_geometry",
  "status": "completed",
  "run_dir": "%s",
  "outputs": {"html_report": "%s"}
}
"""
        % (run, report),
        encoding="utf-8",
    )
    workflow = {item.workflow_id: item for item in default_workflows()}["semiexp_geometry"]

    lines = workflow_state_lines(workflow, tmp_path)

    assert any("latest manifest: completed" in line for line in lines)
    assert any("outputs present: html_report" in line for line in lines)


@pytest.mark.usefixtures("qtbot")
def test_dashboard_lists_workflows(tmp_path, qtbot):
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Merlino 4.0"
    assert _workflow_tree_count(window) == len(default_workflows())
    listed = _workflow_tree_ids(window)
    assert "dvr" in listed
    assert "vpt2_vci" in listed
    assert "semiexp_geometry" in listed
    assert window.menuBar().actions()
    assert {window.menuBar().actions()[i].text() for i in range(window.menuBar().actions().__len__())} >= {
        "Project",
        "Structure",
        "Coordinates",
        "Vibrations",
        "Dynamics",
    }
    window.select_workflow("semiexp_geometry")
    assert window.backend_selector.isEnabled()
    assert not window.semiexp_panel.isHidden()
    assert window.backend_selector.findText("python") >= 0
    assert window.backend_selector.findText("fortran77") >= 0
    assert window.semiexp_preview_table.columnCount() == 5
    window.backend_selector.setCurrentText("fortran77")
    assert window.selected_backends["semiexp_geometry"] == "fortran77"
    assert "selected: fortran77" in window.detail_view.toPlainText()
    window.semiexp_xyz.setText(str(tmp_path / "parent.xyz"))
    window.semiexp_observations.setText(str(tmp_path / "isotopologues.toml"))
    window.semiexp_outdir.setText(str(tmp_path / "semiexp"))
    window.semiexp_fixed.setText("GIC001")
    window.semiexp_qm.setText("GIC002:1.0:0.1:qm")
    window.semiexp_classes.setText("CH:shared:bond(1,2)|bond(1,3);XYH:fixed:angle")
    args = window.semiexp_command_args()
    assert "--backend" in args
    assert "fortran77" in args
    assert args[args.index("--prune-condition") + 1] == "200"
    assert args.count("--parameter-class") == 2
    assert window.semiexp_run_button.isEnabled()
    window.semiexp_iso_table.item(0, 2).setText("1000.0")
    window.semiexp_iso_table.item(0, 3).setText("800.0")
    window.semiexp_iso_table.item(0, 4).setText("600.0")
    toml = window.save_semiexp_observations_toml()
    assert toml.exists()
    assert "A_MHz = 1000.0" in toml.read_text(encoding="utf-8")
    preset = window.save_semiexp_preset()
    assert preset.exists()
    window.semiexp_fixed.clear()
    window.load_semiexp_preset(preset)
    assert window.semiexp_fixed.text() == "GIC001"


@pytest.mark.usefixtures("qtbot")
def test_semiexp_dashboard_preview_validate_and_conditioning(tmp_path, qtbot):
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
    obs = tmp_path / "isotopologues.toml"
    obs.write_text(
        "\n".join(
            [
                "[[isotopologues]]",
                'label = "parent"',
                "[isotopologues.constants]",
                "A_MHz = 822180.189172425",
                "B_MHz = 437776.592728178",
                "C_MHz = 285669.514220614",
                "[isotopologues.vibrational_correction]",
                "delta_A_MHz = 0.0",
                "delta_B_MHz = 0.0",
                "delta_C_MHz = 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)
    window.select_workflow("semiexp_geometry")
    window.semiexp_xyz.setText(str(xyz))
    window.semiexp_observations.setText(str(obs))
    window.semiexp_outdir.setText(str(tmp_path / "run"))

    window.preview_semiexp_gics()
    window.validate_semiexp_input()
    window.preview_semiexp_conditioning()

    text = window.semiexp_command.toPlainText()
    assert window.semiexp_preview_table.rowCount() > 0
    assert "input validation: OK" in text
    assert "condition number" in text


def test_dashboard_launcher_parser_accepts_workdir(tmp_path):
    args = build_parser().parse_args(["--workdir", str(tmp_path)])
    assert args.workdir == tmp_path


def _workflow_tree_count(window):
    total = 0
    for top_idx in range(window.workflow_list.topLevelItemCount()):
        total += window.workflow_list.topLevelItem(top_idx).childCount()
    return total


def _workflow_tree_ids(window):
    ids = set()
    for top_idx in range(window.workflow_list.topLevelItemCount()):
        parent = window.workflow_list.topLevelItem(top_idx)
        for child_idx in range(parent.childCount()):
            ids.add(parent.child(child_idx).data(0, Qt.UserRole))
    return ids
