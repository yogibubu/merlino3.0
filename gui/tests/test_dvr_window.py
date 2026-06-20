import os

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


@pytest.mark.usefixtures("qtbot")
def test_dvr_window_selects_latest_gaussian_log(tmp_path, qtbot):
    old_log = tmp_path / "gauin.log"
    new_log = tmp_path / "scan.out"
    old_log.write_text("old\n", encoding="utf-8")
    new_log.write_text("new\n", encoding="utf-8")
    os.utime(old_log, (1, 1))
    os.utime(new_log, (2, 2))

    window = DVRWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window._select_latest_log(show_message=False)

    assert window.log_edit.text() == str(new_log)
    assert "scan.out" in window.preflight_label.text()


@pytest.mark.usefixtures("qtbot")
def test_dvr_window_preflight_summarizes_gaussian_log(tmp_path, qtbot):
    log = tmp_path / "gauout.log"
    log.write_text(
        "\n".join(
            [
                " Standard orientation:",
                " SCF Done:  E(RB3LYP) = -230.0",
                " QPck001 something",
                " Normal termination of Gaussian 16",
            ]
        ),
        encoding="utf-8",
    )
    window = DVRWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.log_edit.setText(str(log))

    window.preview_log(show_message=False)

    assert "OK:" in window.preflight_label.text()
    assert "SCF energies: 1" in window.output_text.toPlainText()
    assert "puckering GIC markers: 1" in window.output_text.toPlainText()


@pytest.mark.usefixtures("qtbot")
def test_dvr_window_refreshes_result_preview(tmp_path, qtbot):
    outdir = tmp_path / "puckering_dvr_outputs"
    figdir = tmp_path / "puckering_dvr_figs"
    outdir.mkdir()
    figdir.mkdir()
    (outdir / "puckering_dvr_summary.txt").write_text("Lowest levels:\n  0 0.0\n", encoding="utf-8")
    (outdir / "puckering_dvr_levels.csv").write_text("state,energy_cm-1\n0,0.0\n", encoding="utf-8")
    (figdir / "puckering_dvr_potential_levels.pdf").write_bytes(b"%PDF\n")

    window = DVRWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.refresh_results()

    text = window.output_text.toPlainText()
    assert "puckering_dvr_summary.txt" in text
    assert "state,energy_cm-1" in text
    assert "puckering_dvr_potential_levels.pdf" in text
