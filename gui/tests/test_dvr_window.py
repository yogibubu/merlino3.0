import pytest

from advanced.dvr_window import DVRWindow


@pytest.mark.usefixtures("qtbot")
def test_dvr_window_defaults(tmp_path, qtbot):
    window = DVRWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)

    assert "gauin.log" in window.log_edit.text()
    assert window.boundary_combo.currentText() == "periodic"
    assert window.solver_combo.currentText() == "fourier"
    assert window.rotconst_check.isChecked()
    assert not window.cremer_check.isChecked()

