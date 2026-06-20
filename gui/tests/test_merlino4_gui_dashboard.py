from __future__ import annotations

import pytest
from PySide6.QtCore import Qt

from merlino_gui import DashboardWindow, default_workflows
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
    assert "Inputs:" in text
    assert "Outputs:" in text
    assert workflow.service in text


@pytest.mark.usefixtures("qtbot")
def test_dashboard_lists_workflows(tmp_path, qtbot):
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Merlino 4.0"
    assert window.workflow_list.count() == len(default_workflows())
    listed = {
        window.workflow_list.item(i).data(Qt.UserRole)
        for i in range(window.workflow_list.count())
    }
    assert "dvr" in listed
    assert "vpt2_vci" in listed
    assert "semiexp_geometry" in listed
