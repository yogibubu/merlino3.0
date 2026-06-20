from __future__ import annotations

import pytest
from PySide6.QtCore import Qt

from merlino_gui import DashboardWindow, default_workflows
from merlino_gui.app import build_parser
from merlino_gui.dashboard import workflow_detail_text


def test_default_workflows_include_new_scientific_areas():
    workflows = {item.workflow_id: item for item in default_workflows()}
    assert "vpt2_vci" in workflows
    assert "semiexp_geometry" in workflows
    assert workflows["vpt2_vci"].service == "merlino_vpt2_vci"
    assert workflows["semiexp_geometry"].service == "merlino_semiexp"


def test_workflow_detail_text_lists_contract_fields():
    workflow = default_workflows()[0]
    text = workflow_detail_text(workflow)
    assert workflow.title in text
    assert "Category:" in text
    assert "Inputs:" in text
    assert "Outputs:" in text
    assert workflow.service in text


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
