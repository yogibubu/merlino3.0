import pytest
from PySide6.QtWidgets import QLabel

from gui.status_reporter import StatusReporter


@pytest.mark.usefixtures("qtbot")
def test_status_reporter_updates_label_and_summary(tmp_path, qtbot):
    working_dir = tmp_path
    (working_dir / "dos_vib.dat").write_text("x")
    (working_dir / "gui.log").write_text("log")

    label = QLabel()
    qtbot.addWidget(label)

    reporter = StatusReporter(working_dir)
    reporter.update_status_label(
        label=label,
        input_type="harmonic",
        source="xyz",
        vib_q=1.234,
        rovib_q=5.678,
        dos_emin=0.0,
        dos_emax=8000.0,
        dos_bin=50.0,
        dos_T=298.15,
        error=None,
    )

    text = label.text()
    assert "Input: harmonic" in text
    assert "Source: xyz" in text
    assert "Q_vib:" in text
    assert "Q_rovib:" in text
    assert "Emin/Emax/bin:" in text

    reporter.write_summary(
        input_type="harmonic",
        source="xyz",
        dos_emin=0.0,
        dos_emax=8000.0,
        dos_bin=50.0,
        dos_vmax=6,
        dos_ncap="10",
        dos_T=298.15,
        vib_q=1.234,
        rovib_q=5.678,
        error=None,
    )

    summary = (working_dir / "summary.txt").read_text()
    assert "MERLINO GUI SUMMARY" in summary
    assert "Input type: harmonic" in summary
    assert "Source: xyz" in summary
    assert "Q_vib:" in summary
    assert "Q_rovib:" in summary
    assert "dos_vib.dat" in summary
    assert "gui.log" in summary
