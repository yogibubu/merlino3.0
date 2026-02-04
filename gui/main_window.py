from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys


from PySide6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QMessageBox,
    QApplication,
    QPushButton,
    QToolBar,
    QDialog,
    QFormLayout,
    QSpinBox,
    QDoubleSpinBox,
    QLineEdit,
    QDialogButtonBox,
    QCheckBox,
    QLabel,
    QHBoxLayout,
    QTextEdit,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
)
from PySide6.QtCore import Qt, QUrl, QTimer
from PySide6.QtGui import QKeySequence, QAction, QDesktopServices

from .input_panel import InputPanel
from .viewer_panel import ViewerPanel
from .workflow_feedback import ask_open_report
from .project_manager import ProjectManager
from advanced.advanced_window import AdvancedWindow
from .dos_controller import DosController
from .xyzin_service import XyzinService
from .gui_settings import GuiSettings
from .status_reporter import StatusReporter
from .bdpcs3_workflow import run_bdpcs3_report
from geometry.thermo_trasl import read_xyz_from_xyzin
from geometry.vib_anh import read_rotational_block
from .logging_utils import get_gui_logger
from .symmetry_panel import (
    load_symmetry_section,
    parse_equivalent_parameter_classes,
    parse_symmetry_overview,
)
from .similarity_window import SimilarityWindow
from .fragment_pipeline_window import FragmentPipelineWindow


def get_project_root():
    return Path(__file__).resolve().parents[1]


GUI_KEYS = {
    "CHARGE",
    "SPIN_MULTIPLICITY",
    "POINT_GROUP",
    "REPRESENTATION",
    "T_K",
    "P_ATM",
}


class MainWindow(QMainWindow):
    def __init__(self, working_dir: Path = None, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Merlino 3.0")
        self.resize(900, 600)

        self.project_root = get_project_root()
        self.working_dir = Path(working_dir) if working_dir else self.project_root / "working"
        self.working_dir.mkdir(exist_ok=True)

        self.xyzin_path = self.working_dir / "xyzin"
        self._backup_xyzin = self.working_dir / ".xyzin_backup"
        self.advanced_window = None
        self.bdpcs3_version = "legacy"
        self._last_xyzin_stamp = None
        self._refresh_timer = None

        # ✅ Project manager (CREATED, NOT RUN)
        self.manager = ProjectManager(self.working_dir)
        self._input_source = "unknown"
        self._log_path = self.working_dir / "gui.log"
        self.dos_controller = DosController(self.working_dir)
        self.xyzin_service = XyzinService()
        self.status_reporter = StatusReporter(self.working_dir)
        self.logger = get_gui_logger(self._log_path)

        # BASIC defaults
        self.basic_charge = 0
        self.basic_mult = 1
        self.basic_pg = "C1"
        self.basic_T = 298.15
        self.basic_P_atm = 1.0
        self.basic_representation = "Ir"

        # DOS/Q(T) defaults
        self.dos_do_vib = True
        self.dos_do_rovib = True
        self.dos_emin = 0.0
        self.dos_emax = 8000.0
        self.dos_bin = 50.0
        self.dos_vmax = 6
        self.dos_ncap = "10"
        self.dos_T = self.basic_T
        self.dos_emax_rot = None
        self.dos_jmax = None
        self._settings_path = self.working_dir / "gui_settings.json"
        self.gui_settings = None
        self._load_gui_settings()

        self.toolbar = QToolBar("Main")
        self.addToolBar(self.toolbar)

        act_edit_basic = QAction("Edit BASIC", self)
        act_edit_basic.triggered.connect(self.open_basic_dialog)
        self.toolbar.addAction(act_edit_basic)

        self.toolbar.addSeparator()

        act_edit_dos = QAction("DOS/Q(T) Settings", self)
        act_edit_dos.triggered.connect(self.open_dos_dialog)
        self.toolbar.addAction(act_edit_dos)

        self.toolbar.addSeparator()

        act_open_working = QAction("Open working folder", self)
        act_open_working.triggered.connect(self._open_working_folder)
        self.toolbar.addAction(act_open_working)

        self.toolbar.addSeparator()

        act_symmetry = QAction("Symmetry panel", self)
        act_symmetry.triggered.connect(self._open_symmetry_panel)
        self.toolbar.addAction(act_symmetry)
        self.toolbar.addSeparator()
        act_similarity = QAction("Synthon similarity", self)
        act_similarity.triggered.connect(self._open_similarity_window)
        self.toolbar.addAction(act_similarity)
        self.toolbar.addSeparator()
        act_fragment = QAction("Fragment pipeline", self)
        act_fragment.triggered.connect(self._open_fragment_pipeline_window)
        self.toolbar.addAction(act_fragment)

        central = QWidget(self)
        layout = QVBoxLayout(central)

        self.input_panel = InputPanel(
            working_dir=self.working_dir,
            on_update=self.on_xyzin_updated,
        )
        self.viewer_panel = ViewerPanel(self.working_dir)
        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.timeout.connect(self.viewer_panel.refresh)

        self.status_label = QLabel("Status: idle")
        self.status_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.status_label.setTextFormat(Qt.RichText)
        self.status_label.setOpenExternalLinks(False)
        self.status_label.linkActivated.connect(self._open_status_link)
        self.status_label.setToolTip(str(self.working_dir))
        self.status_label.setStyleSheet(
            "color: #4a4a4a; font-size: 11px; padding: 6px 4px;"
        )

        self.banner_label = QLabel("")
        self.banner_label.setStyleSheet(
            "color: #8a4b08; background: #fff3cd; border: 1px solid #ffeeba; "
            "padding: 6px 8px; font-size: 11px;"
        )
        self.banner_label.setVisible(False)

        self.summary_label = QLabel("")
        self.summary_label.setStyleSheet(
            "color: #2f2f2f; font-size: 11px; padding: 4px 6px;"
        )
        self.summary_label.setTextInteractionFlags(Qt.TextSelectableByMouse)

        primary_row = QHBoxLayout()
        primary_row.setSpacing(8)

        primary_basic = QPushButton("Edit BASIC")
        primary_working = QPushButton("Open working")
        primary_similarity = QPushButton("Similarity")
        primary_fragment = QPushButton("Fragment pipeline")

        primary_basic.clicked.connect(self.open_basic_dialog)
        primary_working.clicked.connect(self._open_working_folder)
        primary_similarity.clicked.connect(self._open_similarity_window)
        primary_fragment.clicked.connect(self._open_fragment_pipeline_window)

        primary_basic.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_working.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_similarity.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_fragment.setStyleSheet("font-weight: 600; padding: 6px 12px;")

        primary_row.addWidget(primary_basic)
        primary_row.addWidget(primary_working)
        primary_row.addWidget(primary_similarity)
        primary_row.addWidget(primary_fragment)
        primary_row.addStretch()

        layout.addWidget(self.input_panel)
        layout.addLayout(primary_row)
        layout.addWidget(self.viewer_panel)
        layout.addWidget(self.status_label)
        layout.addWidget(self.banner_label)
        layout.addWidget(self.summary_label)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

        self.setCentralWidget(central)
        self._schedule_refresh()

    # --------------------------------------------------
    # Drag & Drop / Paste
    # --------------------------------------------------
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        if not urls:
            return
        path = Path(urls[0].toLocalFile())
        if path.exists():
            self.input_panel.file_edit.setText(str(path))
            self.input_panel._on_file_commit()

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.Paste):
            text = QApplication.clipboard().text().strip()
            if text:
                path = Path(text)
                if path.exists():
                    self.input_panel.file_edit.setText(str(path))
                    self.input_panel._on_file_commit()
                    return
        super().keyPressEvent(event)

    # --------------------------------------------------
    # XYzin update
    # --------------------------------------------------
    def on_xyzin_updated(self, source="input"):
        self._input_source = source
        self._log_event(f"xyzin updated: source={source}")
        if self.xyzin_path.exists():
            shutil.copy(self.xyzin_path, self._backup_xyzin)

        basic_values = {
            "basic_charge": self.basic_charge,
            "basic_mult": self.basic_mult,
            "basic_pg": self.basic_pg,
            "basic_representation": self.basic_representation,
            "basic_T": self.basic_T,
            "basic_P_atm": self.basic_P_atm,
        }
        self.xyzin_service.update_post_input(
            self.xyzin_path,
            GUI_KEYS,
            basic_values,
            remove_smiles=(source != "smiles"),
        )
        summary = self.xyzin_service.summary(self.xyzin_path)

        reply = QMessageBox.question(
            self,
            "Input updated",
            f"The molecular input has been updated.\n\n"
            f"Source: {source}\n"
            f"{summary}\n\n"
            f"Do you want to continue with the new structure?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )

        if reply == QMessageBox.No:
            if self._backup_xyzin.exists():
                shutil.copy(self._backup_xyzin, self.xyzin_path)
            return

        self._schedule_refresh()

        stamp = self._xyzin_stamp()
        if self._last_xyzin_stamp is not None and stamp == self._last_xyzin_stamp:
            self._log_event("xyzin unchanged; pipelines skipped")
            return

        try:
            self._set_busy(True, "Running workflows…")
            # ✅ CANONICAL WORKFLOW
            self.manager.run_full_workflow()
            self._log_event("pipelines: rotational, thermo, topology (bdpcs3 optional next)")
        except Exception as e:
            QMessageBox.critical(self, "Analysis error", str(e))
            self._log_event(f"analysis error: {e}")
            return
        finally:
            self._set_busy(False)

        try:
            self._run_dos_workflow()
            self._log_event("dos workflow completed")
        except Exception as e:
            QMessageBox.warning(self, "DOS/Q(T) warning", str(e))
            self._update_status(error=str(e))
            self._log_event(f"dos warning: {e}")
        self._last_xyzin_stamp = stamp

        ask_open_report(
            self,
            "Rotational",
            "Rotational OK.",
            "rotational.report",
            working_dir=self.working_dir,
        )
        vib_report = self.working_dir / "vibrational.report"
        if vib_report.exists():
            ask_open_report(
                self,
                "Vibrational",
                "Vibrational OK.",
                "vibrational.report",
                working_dir=self.working_dir,
            )
        ask_open_report(
            self,
            "Thermo",
            "Thermo OK.",
            "thermo.report",
            working_dir=self.working_dir,
        )
        ask_open_report(
            self,
            "Topology",
            "Topology OK. Next step: optional DPCS3 \u2192 BDPCS3.",
            "topology.report",
            working_dir=self.working_dir,
        )
        self._ask_open_symmetry_panel()
        self._ask_run_bdpcs3_after_topology()
        self._update_fchkin_banner()
        self._ask_open_advanced()

    # --------------------------------------------------
    # BASIC dialog
    # --------------------------------------------------
    def open_basic_dialog(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Edit BASIC")
        layout = QFormLayout(dlg)

        spin_charge = QSpinBox(dlg)
        spin_charge.setRange(-10, 10)
        spin_charge.setValue(self.basic_charge)

        spin_mult = QSpinBox(dlg)
        spin_mult.setRange(1, 20)
        spin_mult.setValue(self.basic_mult)

        edit_pg = QLineEdit(self.basic_pg)

        spin_T = QDoubleSpinBox(dlg)
        spin_T.setRange(0.0, 5000.0)
        spin_T.setDecimals(2)
        spin_T.setValue(self.basic_T)

        spin_P = QDoubleSpinBox(dlg)
        spin_P.setRange(0.0, 100000.0)
        spin_P.setDecimals(6)
        spin_P.setValue(self.basic_P_atm)

        edit_rep = QLineEdit(self.basic_representation)

        layout.addRow("Charge", spin_charge)
        layout.addRow("Multiplicity", spin_mult)
        layout.addRow("Point group", edit_pg)
        layout.addRow("Temperature (K)", spin_T)
        layout.addRow("Pressure (atm)", spin_P)
        layout.addRow("Representation", edit_rep)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addRow(buttons)

        reset_btn = QPushButton("Reset to defaults")
        def _reset_defaults():
            chk_vib.setChecked(True)
            chk_rovib.setChecked(True)
            spin_emin.setValue(0.0)
            spin_emax.setValue(8000.0)
            spin_bin.setValue(50.0)
            spin_vmax.setValue(6)
            edit_ncap.setText("10")
            spin_T.setValue(self.basic_T)
            spin_emax_rot.setValue(0.0)
            spin_jmax.setValue(0)
        reset_btn.clicked.connect(_reset_defaults)
        layout.addRow(reset_btn)

        if dlg.exec() != QDialog.Accepted:
            return

        self.basic_charge = spin_charge.value()
        self.basic_mult = spin_mult.value()
        self.basic_pg = edit_pg.text().strip() or "C1"
        self.basic_T = spin_T.value()
        self.basic_P_atm = spin_P.value()

        rep = edit_rep.text().strip() or "Ir"
        if rep not in {"Ir", "IIr", "IIIr", "Il", "IIl", "IIIl"}:
            QMessageBox.warning(self, "Invalid Representation", "Reset to Ir.")
            rep = "Ir"
        self.basic_representation = rep

        self._update_basic_section_in_xyzin()
        self._schedule_refresh()

    # --------------------------------------------------
    # DOS/Q(T) dialog
    # --------------------------------------------------
    def open_dos_dialog(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("DOS/Q(T) Settings")
        layout = QFormLayout(dlg)

        chk_vib = QCheckBox("Compute vib DOS/Q(T)")
        chk_vib.setChecked(self.dos_do_vib)

        chk_rovib = QCheckBox("Compute rovib DOS/Q(T)")
        chk_rovib.setChecked(self.dos_do_rovib)
        if not self._has_rotational_block():
            chk_rovib.setChecked(False)
            chk_rovib.setEnabled(False)
            chk_rovib.setToolTip("Missing #ROTATIONAL: rovib not available")

        spin_emin = QDoubleSpinBox(dlg)
        spin_emin.setRange(0.0, 1.0e6)
        spin_emin.setDecimals(2)
        spin_emin.setValue(self.dos_emin)

        spin_emax = QDoubleSpinBox(dlg)
        spin_emax.setRange(0.0, 1.0e6)
        spin_emax.setDecimals(2)
        spin_emax.setValue(self.dos_emax)

        spin_bin = QDoubleSpinBox(dlg)
        spin_bin.setRange(0.1, 1.0e6)
        spin_bin.setDecimals(2)
        spin_bin.setValue(self.dos_bin)

        spin_vmax = QSpinBox(dlg)
        spin_vmax.setRange(0, 200)
        spin_vmax.setValue(self.dos_vmax)

        edit_ncap = QLineEdit(self.dos_ncap)

        spin_T = QDoubleSpinBox(dlg)
        spin_T.setRange(0.1, 5000.0)
        spin_T.setDecimals(2)
        spin_T.setValue(self.dos_T)

        spin_emax_rot = QDoubleSpinBox(dlg)
        spin_emax_rot.setRange(0.0, 1.0e6)
        spin_emax_rot.setDecimals(2)
        spin_emax_rot.setValue(self.dos_emax_rot or 0.0)

        spin_jmax = QSpinBox(dlg)
        spin_jmax.setRange(0, 5000)
        spin_jmax.setValue(self.dos_jmax or 0)

        layout.addRow(chk_vib)
        layout.addRow(chk_rovib)
        layout.addRow("Emin (cm-1)", spin_emin)
        layout.addRow("Emax (cm-1)", spin_emax)
        layout.addRow("Bin (cm-1)", spin_bin)
        layout.addRow("vmax", spin_vmax)
        layout.addRow("ncap (comma or scalar)", edit_ncap)
        layout.addRow("Temperature (K)", spin_T)
        layout.addRow("Emax rot (cm-1, 0=auto)", spin_emax_rot)
        layout.addRow("Jmax (0=auto)", spin_jmax)

        err_label = QLabel("")
        err_label.setStyleSheet("color: #b00020;")
        layout.addRow(err_label)

        def _validate_live():
            msg = None
            if spin_emax.value() <= spin_emin.value():
                msg = "Emax must be > Emin."
            elif spin_bin.value() <= 0.0:
                msg = "Bin must be > 0."
            elif spin_T.value() <= 0.0:
                msg = "T must be > 0."
            elif spin_vmax.value() <= 0:
                msg = "vmax must be > 0."
            err_label.setText(msg or "")
            return msg is None

        spin_emin.valueChanged.connect(_validate_live)
        spin_emax.valueChanged.connect(_validate_live)
        spin_bin.valueChanged.connect(_validate_live)
        spin_T.valueChanged.connect(_validate_live)
        spin_vmax.valueChanged.connect(_validate_live)
        _validate_live()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addRow(buttons)

        if dlg.exec() != QDialog.Accepted:
            return

        if not _validate_live():
            QMessageBox.warning(self, "Invalid DOS settings", err_label.text())
            return
        if spin_emax.value() <= spin_emin.value():
            QMessageBox.warning(self, "Invalid DOS range", "Emax must be > Emin.")
            return
        if spin_bin.value() <= 0.0:
            QMessageBox.warning(self, "Invalid bin", "Bin must be > 0.")
            return
        if spin_T.value() <= 0.0:
            QMessageBox.warning(self, "Invalid temperature", "T must be > 0.")
            return
        if spin_vmax.value() <= 0:
            QMessageBox.warning(self, "Invalid vmax", "vmax must be > 0.")
            return

        self.dos_do_vib = chk_vib.isChecked()
        self.dos_do_rovib = chk_rovib.isChecked()
        self.dos_emin = spin_emin.value()
        self.dos_emax = spin_emax.value()
        self.dos_bin = spin_bin.value()
        self.dos_vmax = spin_vmax.value()
        self.dos_ncap = edit_ncap.text().strip()
        self.dos_T = spin_T.value()
        self.dos_emax_rot = spin_emax_rot.value() or None
        self.dos_jmax = spin_jmax.value() or None
        self._save_gui_settings()

    def _update_basic_section_in_xyzin(self):
        basic_values = {
            "basic_charge": self.basic_charge,
            "basic_mult": self.basic_mult,
            "basic_pg": self.basic_pg,
            "basic_representation": self.basic_representation,
            "basic_T": self.basic_T,
            "basic_P_atm": self.basic_P_atm,
        }
        self.xyzin_service.update_basic_section(self.xyzin_path, GUI_KEYS, basic_values)

    def _ask_open_advanced(self):
        reply = QMessageBox.question(
            self,
            "Preparation completed",
            "Preparation completed successfully.\n\n"
            "Do you want to continue with Advanced calculations?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )

        if reply != QMessageBox.Yes:
            return

        if self.advanced_window is None:
            self.advanced_window = AdvancedWindow(self.working_dir, parent=self)

        self.advanced_window.show()
        self.advanced_window.raise_()
        self.advanced_window.activateWindow()

    # --------------------------------------------------
    # DOS/Q(T) workflow
    # --------------------------------------------------
    def _run_dos_workflow(self):
        settings = {
            "dos_do_vib": self.dos_do_vib,
            "dos_do_rovib": self.dos_do_rovib,
            "dos_emin": self.dos_emin,
            "dos_emax": self.dos_emax,
            "dos_bin": self.dos_bin,
            "dos_vmax": self.dos_vmax,
            "dos_ncap": self.dos_ncap,
            "dos_T": self.dos_T,
            "dos_emax_rot": self.dos_emax_rot,
            "dos_jmax": self.dos_jmax,
        }
        result = self.dos_controller.run(self.xyzin_path, settings)
        if result is None:
            return
        self._update_status(
            input_type=result.get("input_type"),
            vib_q=result.get("vib_q"),
            rovib_q=result.get("rovib_q"),
            error=result.get("error"),
        )

    def _update_status(self, input_type=None, vib_q=None, rovib_q=None, error=None):
        self.status_reporter.update_status_label(
            self.status_label,
            input_type,
            self._input_source,
            vib_q,
            rovib_q,
            self.dos_emin,
            self.dos_emax,
            self.dos_bin,
            self.dos_T,
            error,
        )
        self._update_results_panel(vib_q, rovib_q)
        self.status_reporter.write_summary(
            input_type,
            self._input_source,
            self.dos_emin,
            self.dos_emax,
            self.dos_bin,
            self.dos_vmax,
            self.dos_ncap,
            self.dos_T,
            vib_q,
            rovib_q,
            error,
        )
        if error:
            self.banner_label.setText(error)
            self.banner_label.setVisible(True)
        else:
            if not self.banner_label.text():
                self.banner_label.setVisible(False)

    def _update_results_panel(self, vib_q=None, rovib_q=None):
        parts = []
        if vib_q is not None:
            parts.append(f"Q_vib: {vib_q:.6e}")
        if rovib_q is not None:
            parts.append(f"Q_rovib: {rovib_q:.6e}")
        rot = {}
        if self.xyzin_path.exists():
            try:
                rot = read_rotational_block(str(self.xyzin_path))
            except Exception:
                rot = {}
        A = rot.get("a_mhz")
        B = rot.get("b_mhz")
        C = rot.get("c_mhz")
        if A and B and C:
            parts.append(f"A/B/C (MHz): {A} / {B} / {C}")
        self.summary_label.setText(" | ".join(parts))

    def _set_busy(self, busy: bool, text: str | None = None):
        self.toolbar.setEnabled(not busy)
        self.input_panel.setEnabled(not busy)
        if busy:
            self.status_label.setText(text or "Running…")
        else:
            if not self.status_label.text():
                self.status_label.setText("Status: idle")

    def _schedule_refresh(self, delay_ms: int = 120):
        if self._refresh_timer is None:
            return
        self._refresh_timer.start(delay_ms)

    def _xyzin_stamp(self):
        if not self.xyzin_path.exists():
            return None
        try:
            stat = self.xyzin_path.stat()
        except OSError:
            return None
        return (stat.st_mtime, stat.st_size)

    def _update_fchkin_banner(self):
        fchkin = self.working_dir / "fchkin"
        if not fchkin.exists() or not self.xyzin_path.exists():
            self.banner_label.setVisible(False)
            return
        try:
            from geometry.vibrational import read_fchk_masses
            symbols, _coords = read_xyz_from_xyzin(str(self.xyzin_path))
            masses = read_fchk_masses(str(fchkin))
            if len(masses) != len(symbols):
                msg = (
                    "Vibrations skipped: fchkin atom count does not match xyzin "
                    f"({len(masses)} vs {len(symbols)})"
                )
                self.banner_label.setText(msg)
                self.banner_label.setVisible(True)
            else:
                self.banner_label.setVisible(False)
        except Exception:
            return

    def _load_gui_settings(self):
        defaults = GuiSettings(
            dos_do_vib=self.dos_do_vib,
            dos_do_rovib=self.dos_do_rovib,
            dos_emin=self.dos_emin,
            dos_emax=self.dos_emax,
            dos_bin=self.dos_bin,
            dos_vmax=self.dos_vmax,
            dos_ncap=self.dos_ncap,
            dos_T=self.dos_T,
            dos_emax_rot=self.dos_emax_rot,
            dos_jmax=self.dos_jmax,
        )
        self.gui_settings = GuiSettings.from_json(self._settings_path, defaults=defaults)
        self.dos_do_vib = bool(self.gui_settings.dos_do_vib)
        self.dos_do_rovib = bool(self.gui_settings.dos_do_rovib)
        self.dos_emin = float(self.gui_settings.dos_emin)
        self.dos_emax = float(self.gui_settings.dos_emax)
        self.dos_bin = float(self.gui_settings.dos_bin)
        self.dos_vmax = int(self.gui_settings.dos_vmax)
        self.dos_ncap = str(self.gui_settings.dos_ncap)
        self.dos_T = float(self.gui_settings.dos_T)
        self.dos_emax_rot = self.gui_settings.dos_emax_rot
        self.dos_jmax = self.gui_settings.dos_jmax

    def _save_gui_settings(self):
        self.gui_settings = GuiSettings(
            dos_do_vib=self.dos_do_vib,
            dos_do_rovib=self.dos_do_rovib,
            dos_emin=self.dos_emin,
            dos_emax=self.dos_emax,
            dos_bin=self.dos_bin,
            dos_vmax=self.dos_vmax,
            dos_ncap=self.dos_ncap,
            dos_T=self.dos_T,
            dos_emax_rot=self.dos_emax_rot,
            dos_jmax=self.dos_jmax,
        )
        try:
            self.gui_settings.to_json(self._settings_path)
        except Exception:
            pass

    def _open_working_folder(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.working_dir)))

    def _open_text_file(self, path: Path):
        path = Path(path)
        if not path.exists():
            QMessageBox.warning(self, "File not found", str(path))
            return
        try:
            if sys.platform == "darwin":
                subprocess.Popen(["open", "-a", "TextEdit", str(path)])
            elif sys.platform.startswith("win"):
                subprocess.Popen(["notepad.exe", str(path)])
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except Exception:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def _ask_open_symmetry_panel(self):
        report_path = self.working_dir / "topology.report"
        if not report_path.exists():
            return
        reply = QMessageBox.question(
            self,
            "Symmetry Analysis",
            "Topology includes point group and equivalent-parameter classes.\n\n"
            "Do you want to open the dedicated symmetry panel?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )
        if reply == QMessageBox.Yes:
            self._open_symmetry_panel()

    def _open_similarity_window(self):
        dlg = SimilarityWindow(self.working_dir, self)
        dlg.exec()

    def _open_fragment_pipeline_window(self):
        dlg = FragmentPipelineWindow(self.working_dir, self)
        dlg.exec()

    def _open_symmetry_panel(self):
        report_path = self.working_dir / "topology.report"
        sym_text = load_symmetry_section(report_path)
        if not sym_text:
            QMessageBox.information(
                self,
                "Symmetry panel",
                "No symmetry section found in topology.report.\n"
                "Run topology analysis first.",
            )
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Symmetry Panel")
        dlg.resize(760, 520)
        layout = QVBoxLayout(dlg)

        title = QLabel(
            "Point group and equivalent internal-parameter classes\n"
            "(from topology.report)"
        )
        title.setStyleSheet("font-weight: 600;")
        layout.addWidget(title)

        tabs = QTabWidget(dlg)
        layout.addWidget(tabs)

        overview = parse_symmetry_overview(sym_text)
        pg = overview.get("point_group") or "n/a"
        atom_classes = overview.get("atom_classes") or []
        overview_text = QTextEdit(dlg)
        overview_text.setReadOnly(True)
        overview_text.setPlainText(
            f"Point group: {pg}\n\n"
            "Equivalent atom classes:\n"
            + ("\n".join(atom_classes) if atom_classes else "none")
        )
        tabs.addTab(overview_text, "Overview")

        classes_by_kind = parse_equivalent_parameter_classes(sym_text)
        for kind in ("bond", "angle", "dihedral", "out_of_plane", "linear_bend"):
            table = QTableWidget(dlg)
            table.setColumnCount(2)
            table.setHorizontalHeaderLabels(["Class", "Equivalent members"])
            rows = classes_by_kind.get(kind, [])
            table.setRowCount(max(1, len(rows)))
            if rows:
                for r, (cid, members) in enumerate(rows):
                    table.setItem(r, 0, QTableWidgetItem(cid))
                    table.setItem(r, 1, QTableWidgetItem(members))
            else:
                table.setItem(0, 0, QTableWidgetItem("-"))
                table.setItem(0, 1, QTableWidgetItem("No equivalent classes"))
            hdr = table.horizontalHeader()
            hdr.setSectionResizeMode(0, QHeaderView.ResizeToContents)
            hdr.setSectionResizeMode(1, QHeaderView.Stretch)
            table.verticalHeader().setVisible(False)
            table.setEditTriggers(QTableWidget.NoEditTriggers)
            table.setSelectionBehavior(QTableWidget.SelectRows)
            tabs.addTab(table, kind)

        raw_text = QTextEdit(dlg)
        raw_text.setReadOnly(True)
        raw_text.setPlainText(sym_text)
        tabs.addTab(raw_text, "Raw report")

        btn_row = QHBoxLayout()
        layout.addLayout(btn_row)

        btn_open_report = QPushButton("Open topology.report")
        btn_copy = QPushButton("Copy panel text")
        btn_close = QPushButton("Close")

        btn_row.addWidget(btn_open_report)
        btn_row.addWidget(btn_copy)
        btn_row.addStretch()
        btn_row.addWidget(btn_close)

        btn_open_report.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(report_path)))
        )
        btn_copy.clicked.connect(lambda: QApplication.clipboard().setText(sym_text))
        btn_close.clicked.connect(dlg.accept)

        dlg.exec()

    def _generate_bdpcs3_report(self):
        result = run_bdpcs3_report(
            self,
            self.working_dir,
            self.xyzin_path,
            self.basic_charge,
            self.basic_mult,
            self.bdpcs3_version,
        )
        if result is None:
            return
        self.bdpcs3_version = result.get("version", self.bdpcs3_version)
        self._log_event("bdpcs3 report generated")
        self._show_bdpcs3_actions(result)
        ask_open_report(
            self,
            "BDPCS3",
            "BDPCS3 report generated.",
            "bdpcs3.report",
            working_dir=self.working_dir,
        )

    def _ask_run_bdpcs3_after_topology(self):
        reply = QMessageBox.question(
            self,
            "Run BDPCS3",
            "Topology analysis completed.\n\n"
            "Do you want to run DPCS3 \u2192 BDPCS3 now?\n"
            "This will also generate the Gaussian input file (BDPCS3.gjf).",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )
        if reply != QMessageBox.Yes:
            return
        self._generate_bdpcs3_report()

    def _has_rotational_block(self):
        return self.xyzin_service.has_rotational_block(self.xyzin_path)

    def _show_bdpcs3_actions(self, result: dict):
        report_path = result.get("report_path")
        gjf_path = result.get("gjf_path")

        dlg = QDialog(self)
        dlg.setWindowTitle("BDPCS3 outputs")
        layout = QVBoxLayout(dlg)

        info = []
        if report_path:
            info.append(f"Report: {report_path}")
        if gjf_path:
            info.append(f"Gaussian input: {gjf_path}")
        label = QLabel("\n".join(info) if info else "Files generated.")
        label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(label)

        btn_row = QHBoxLayout()
        layout.addLayout(btn_row)

        btn_open_report = QPushButton("Open report")
        btn_open_gjf = QPushButton("Open BDPCS3.gjf (text)")
        btn_copy_report = QPushButton("Copy report path")
        btn_copy_gjf = QPushButton("Copy gjf path")

        btn_row.addWidget(btn_open_report)
        btn_row.addWidget(btn_open_gjf)
        btn_row.addWidget(btn_copy_report)
        btn_row.addWidget(btn_copy_gjf)

        btn_open_report.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(report_path)))
        )
        btn_open_gjf.clicked.connect(
            lambda: self._open_text_file(gjf_path)
        )
        btn_copy_report.clicked.connect(
            lambda: QApplication.clipboard().setText(str(report_path))
        )
        btn_copy_gjf.clicked.connect(
            lambda: QApplication.clipboard().setText(str(gjf_path))
        )

        if not report_path:
            btn_open_report.setEnabled(False)
            btn_copy_report.setEnabled(False)
        if not gjf_path:
            btn_open_gjf.setEnabled(False)
            btn_copy_gjf.setEnabled(False)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.accept)
        layout.addWidget(close_btn)

        dlg.exec()

    def _log_event(self, msg: str):
        try:
            self.logger.info(msg)
        except Exception:
            pass

    def _open_status_link(self, url: str):
        if not url:
            return
        QDesktopServices.openUrl(QUrl(url))
