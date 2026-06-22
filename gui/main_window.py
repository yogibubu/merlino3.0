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
    QRadioButton,
    QLabel,
    QHBoxLayout,
    QTextEdit,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QButtonGroup,
)
from PySide6.QtCore import Qt, QUrl, QTimer
from PySide6.QtGui import QKeySequence, QAction, QDesktopServices

from .input_panel import InputPanel
from .viewer_panel import ViewerPanel
from .workflow_feedback import ask_open_report
from .project_manager import ProjectManager
from advanced.advanced_window import AdvancedWindow
from advanced.dvr_window import DVRWindow
from .dos_controller import DosController
from .xyzin_service import XyzinService
from .gui_settings import GuiSettings
from .status_reporter import StatusReporter
from .bdpcs3_workflow import run_bdpcs3_report
from geometry.thermo_trasl import read_xyz_from_xyzin
from geometry.thermo_pipeline import read_xyz_from_xyzin as read_xyz_full_from_xyzin
from geometry.vib_anh import read_rotational_block
from .logging_utils import get_gui_logger
from .symmetry_panel import (
    load_symmetry_section,
    parse_equivalent_parameter_classes,
    parse_symmetry_overview,
)
from .similarity_window import SimilarityWindow
from .fragment_pipeline_window import FragmentPipelineWindow
from topology.elements import atomic_number
from geometry.isotopes_table import get_default_isotope, get_isotopes
from merlino_core.isotopologues import XyzinIsotopologueRecord, merge_xyzin_isotopologue_records
from merlino_fit.survibfit.modify_geom import write_xyz, read_xyz
from merlino_fit.survibfit.fragment_pipeline import run_fragment_pipeline, write_fragment_view_html
from merlino_fit.survibfit.fragment_delta_correction import (
    prepare_hpcs2_delta_workflow,
    apply_delta_correction,
)


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
        self.dvr_window = None
        self.bdpcs3_version = "updated"
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
        self.toolbar.addSeparator()
        act_dvr = QAction("DVR", self)
        act_dvr.triggered.connect(self._open_dvr_window)
        self.toolbar.addAction(act_dvr)
        self.toolbar.addSeparator()
        act_isot = QAction("Isotopologues…", self)
        act_isot.triggered.connect(self._open_isotopologues_dialog)
        self.toolbar.addAction(act_isot)
        self.toolbar.addSeparator()
        act_actions = QAction("Run action…", self)
        act_actions.triggered.connect(self._run_action_dialog)
        self.toolbar.addAction(act_actions)

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
        primary_actions = QPushButton("Run action…")

        primary_basic.clicked.connect(self.open_basic_dialog)
        primary_working.clicked.connect(self._open_working_folder)
        primary_similarity.clicked.connect(self._open_similarity_window)
        primary_fragment.clicked.connect(self._open_fragment_pipeline_window)
        primary_actions.clicked.connect(self._run_action_dialog)

        primary_basic.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_working.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_similarity.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_fragment.setStyleSheet("font-weight: 600; padding: 6px 12px;")
        primary_actions.setStyleSheet("font-weight: 600; padding: 6px 12px;")

        primary_row.addWidget(primary_basic)
        primary_row.addWidget(primary_working)
        primary_row.addWidget(primary_similarity)
        primary_row.addWidget(primary_fragment)
        primary_row.addWidget(primary_actions)
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
        self._startup_hidden = True
        self.hide()
        QTimer.singleShot(0, self._show_startup_dialog)

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
        action = self._prompt_action_choice()
        if not action:
            return
        self._run_selected_action(action)
        self._last_xyzin_stamp = self._xyzin_stamp()

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

    def _open_dvr_window(self):
        if self.dvr_window is None:
            self.dvr_window = DVRWindow(
                self.working_dir,
                get_project_root(),
                parent=self,
            )
        self.dvr_window.show()
        self.dvr_window.raise_()
        self.dvr_window.activateWindow()

    # --------------------------------------------------
    # Action menu
    # --------------------------------------------------
    def _prompt_action_choice(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Select action")
        layout = QVBoxLayout(dlg)

        title = QLabel("Choose what to run for the current input:")
        title.setStyleSheet("font-weight: 600;")
        layout.addWidget(title)

        group = QButtonGroup(dlg)
        group.setExclusive(True)
        options = [
            ("geometry", "Geometry only (no analysis)"),
            ("symmetry", "Symmetry (runs topology, then opens symmetry panel)"),
            ("rotational", "Rotational"),
            ("vibrational", "Vibrational"),
            ("thermo", "Thermo"),
            ("dos", "DOS/Q(T)"),
            ("topology", "Topology"),
            ("bdpcs3", "BDPCS3"),
            ("hpcs2", "HPCS2 (PCS2 geometry from HPCS2 base)"),
        ]
        btns = {}
        for idx, (key, label) in enumerate(options):
            btn = QRadioButton(label, dlg)
            btns[key] = btn
            group.addButton(btn, idx)
            layout.addWidget(btn)

        btns["geometry"].setChecked(True)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        if dlg.exec() != QDialog.Accepted:
            return None

        for key, btn in btns.items():
            if btn.isChecked():
                return key
        return None

    def _export_xyzin_xyz(self) -> Path | None:
        if not self.xyzin_path.exists():
            return None
        try:
            _nat, comment, symbols, coords, _tail = read_xyz_full_from_xyzin(
                str(self.xyzin_path)
            )
        except Exception:
            return None
        out = self.working_dir / f"{self.xyzin_path.stem}.xyz"
        try:
            write_xyz(out, symbols, coords, comment=comment or "xyzin export")
        except Exception:
            return None
        return out

    def _resolve_library_dirs(self):
        se_local = self.working_dir / "projects" / "se_library"
        pcs2_local = self.working_dir / "projects" / "pcs2_library"
        hpcs2_local = self.working_dir / "projects" / "hpcs2_library"
        se_parent = self.working_dir.parent / "projects" / "se_library"
        pcs2_parent = self.working_dir.parent / "projects" / "pcs2_library"
        hpcs2_parent = self.working_dir.parent / "projects" / "hpcs2_library"
        se_dir = se_local if se_local.is_dir() else se_parent
        pcs2_dir = pcs2_local if pcs2_local.is_dir() else pcs2_parent
        hpcs2_dir = hpcs2_local if hpcs2_local.is_dir() else hpcs2_parent
        return se_dir, pcs2_dir, hpcs2_dir

    def _run_hpcs2_flow(self, query_xyz: Path):
        se_dir, pcs2_dir, hpcs2_dir = self._resolve_library_dirs()
        if not pcs2_dir.is_dir():
            QMessageBox.warning(self, "HPCS2", "PCS2 library directory not found.")
            return
        if not hpcs2_dir.is_dir():
            QMessageBox.warning(self, "HPCS2", "HPCS2 library directory not found.")
            return
        symm_tol = self._prompt_symmetry_tolerance()
        if symm_tol is None:
            return
        symm_apply = False
        apply_choice = self._prompt_symmetrize_coords()
        if apply_choice is None:
            return
        symm_apply = bool(apply_choice)
        try:
            self.manager.run_topology(
                symm_tol=float(symm_tol),
                symmetrize_coords=symm_apply,
            )
        except Exception as e:
            QMessageBox.warning(self, "HPCS2", f"Topology failed: {e}")
            return
        ask_open_report(
            self,
            "Topology",
            "Topology OK.",
            "topology.report",
            working_dir=self.working_dir,
        )
        if symm_apply:
            self._apply_symmetrized_xyz()
        out_dir = self.working_dir / "fragment_reports"
        report = run_fragment_pipeline(
            query_xyz,
            se_dir,
            pcs2_dir,
            out_dir,
            use_pcs2_only=True,
            library_label="PCS2/HPCS2",
        )
        view_path = out_dir / "fragment_view.html"
        if view_path.exists():
            reply = QMessageBox.question(
                self,
                "Open fragment viewer",
                "Open 3D fragment viewer now?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes,
            )
            if reply == QMessageBox.Yes:
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(view_path)))
        manifest = prepare_hpcs2_delta_workflow(
            query_xyz,
            out_dir / "fragment_pipeline.json",
            hpcs2_dir,
            out_dir / "delta_bundle_hpcs2",
        )
        if not manifest.get("entries"):
            QMessageBox.warning(
                self,
                "HPCS2",
                "No valid PCS2/HPCS2 fragment pairs found. Check libraries and rerun.",
            )
            return
        out_xyz = self.working_dir / f"{query_xyz.stem}.pcs2.xyz"
        meta = apply_delta_correction(
            query_xyz,
            out_dir / "delta_bundle_hpcs2" / "delta_manifest.json",
            out_xyz,
        )
        write_fragment_view_html(
            query_xyz,
            report.get("fragments", []),
            out_dir,
            corrected_xyz=out_xyz,
            corrected_label="PCS2",
        )
        QMessageBox.information(
            self,
            "HPCS2",
            f"PCS2 geometry written:\n{out_xyz}\n\n"
            f"Fragments used: {meta.get('entries_used', 0)} / {meta.get('entries_total', 0)}",
        )

    def _apply_symmetrized_xyz(self):
        symm_xyz = self.working_dir / "symmetrized.xyz"
        if not symm_xyz.exists():
            return
        try:
            atoms, coords, _ = read_xyz(symm_xyz)
            from .xyzin_utils import replace_xyz_block

            lines = [str(len(atoms)), "Symmetrized coordinates"]
            for a, (x, y, z) in zip(atoms, coords):
                lines.append(f"{a} {x: .6f} {y: .6f} {z: .6f}")
            replace_xyz_block(lines)
            self._log_event("xyzin updated from symmetrized.xyz")
        except Exception:
            return

    def _run_selected_action(self, action: str):
        try:
            self._set_busy(True, f"Running {action}…")
            if action == "geometry":
                self._schedule_refresh()
                self.status_label.setText("Geometry ready.")
                self._log_event("action: geometry only")
                return
            if action == "rotational":
                self.manager.run_rotational()
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
                self._log_event("action: rotational")
                return
            if action == "vibrational":
                self.manager.run_vibrational()
                if (self.working_dir / "vibrational.report").exists():
                    ask_open_report(
                        self,
                        "Vibrational",
                        "Vibrational OK.",
                        "vibrational.report",
                        working_dir=self.working_dir,
                    )
                else:
                    QMessageBox.information(
                        self,
                        "Vibrational",
                        "No vibrational data available (missing fchkin and #VIBRATIONAL block).",
                    )
                self._log_event("action: vibrational")
                return
            if action == "thermo":
                self.manager.run_thermo()
                ask_open_report(
                    self,
                    "Thermo",
                    "Thermo OK.",
                    "thermo.report",
                    working_dir=self.working_dir,
                )
                self._log_event("action: thermo")
                return
            if action == "dos":
                self._run_dos_workflow()
                self._log_event("action: dos")
                return
            if action == "topology":
                symm_tol = self._prompt_symmetry_tolerance()
                if symm_tol is None:
                    return
                apply_choice = self._prompt_symmetrize_coords()
                if apply_choice is None:
                    return
                symm_apply = bool(apply_choice)
                self.manager.run_topology(symm_tol=float(symm_tol), symmetrize_coords=symm_apply)
                if symm_apply:
                    self._apply_symmetrized_xyz()
                ask_open_report(
                    self,
                    "Topology",
                    "Topology OK.",
                    "topology.report",
                    working_dir=self.working_dir,
                )
                self._log_event("action: topology")
                return
            if action == "symmetry":
                symm_tol = self._prompt_symmetry_tolerance()
                if symm_tol is None:
                    return
                apply_choice = self._prompt_symmetrize_coords()
                if apply_choice is None:
                    return
                symm_apply = bool(apply_choice)
                self.manager.run_topology(symm_tol=float(symm_tol), symmetrize_coords=symm_apply)
                if symm_apply:
                    self._apply_symmetrized_xyz()
                self._open_symmetry_panel()
                self._log_event("action: symmetry")
                return
            if action == "bdpcs3":
                self._generate_bdpcs3_report()
                self._log_event("action: bdpcs3")
                return
            if action == "hpcs2":
                query_xyz = self._export_xyzin_xyz()
                if query_xyz is None or not query_xyz.exists():
                    QMessageBox.warning(
                        self,
                        "HPCS2",
                        "Unable to export current xyzin to XYZ for HPCS2 workflow.",
                    )
                    return
                self._run_hpcs2_flow(query_xyz)
                self._log_event("action: hpcs2")
                return
        except Exception as e:
            QMessageBox.critical(self, "Analysis error", str(e))
            self._log_event(f"analysis error: {e}")
        finally:
            self._set_busy(False)
            self._ask_run_another_action()

    def _run_action_dialog(self):
        action = self._prompt_action_choice()
        if not action:
            return
        self._run_selected_action(action)

    def _ask_run_another_action(self):
        reply = QMessageBox.question(
            self,
            "Run another action",
            "Do you want to run another action on the same input?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )
        if reply != QMessageBox.Yes:
            return
        action = self._prompt_action_choice()
        if not action:
            return
        self._run_selected_action(action)

    # --------------------------------------------------
    # Startup flow
    # --------------------------------------------------
    def _show_startup_dialog(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Select input type")
        layout = QVBoxLayout(dlg)
        title = QLabel("Choose input type to begin:")
        title.setStyleSheet("font-weight: 600;")
        layout.addWidget(title)

        group = QButtonGroup(dlg)
        group.setExclusive(True)
        options = ["SMILES", "XYZ", "Gaussian", "Molpro", "MRCC"]
        btns = {}
        for idx, name in enumerate(options):
            btn = QRadioButton(name, dlg)
            btns[name] = btn
            group.addButton(btn, idx)
            layout.addWidget(btn)
        btns["SMILES"].setChecked(True)

        rep_row = QHBoxLayout()
        rep_label = QLabel("Representation:")
        rep_edit = QLineEdit(self.basic_representation)
        rep_edit.setPlaceholderText("Ir / IIr / IIIr / Il / IIl / IIIl")
        rep_row.addWidget(rep_label)
        rep_row.addWidget(rep_edit)
        layout.addLayout(rep_row)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        if dlg.exec() != QDialog.Accepted:
            self.close()
            return

        rep = rep_edit.text().strip() or "Ir"
        if rep not in {"Ir", "IIr", "IIIr", "Il", "IIl", "IIIl"}:
            QMessageBox.warning(self, "Invalid Representation", "Reset to Ir.")
            rep = "Ir"
        self.basic_representation = rep
        self._update_basic_section_in_xyzin()

        selected = "SMILES"
        for name, btn in btns.items():
            if btn.isChecked():
                selected = name
                break
        self.input_panel.set_input_type(selected)
        if selected != "Gaussian":
            self.input_panel.open_file_dialog()
        if self._startup_hidden:
            self._startup_hidden = False
            self.show()

    def _prompt_symmetry_tolerance(self) -> float | None:
        if self._input_source == "smiles":
            return 5.0e-2
        dlg = QDialog(self)
        dlg.setWindowTitle("Symmetry tolerance")
        layout = QVBoxLayout(dlg)
        title = QLabel("Choose symmetry tolerance:")
        title.setStyleSheet("font-weight: 600;")
        layout.addWidget(title)

        group = QButtonGroup(dlg)
        group.setExclusive(True)
        opt_strict = QRadioButton("Strict (1e-3 Å)", dlg)
        opt_loose = QRadioButton("Loose (5e-2 Å)", dlg)
        opt_strict.setChecked(True)
        group.addButton(opt_strict, 0)
        group.addButton(opt_loose, 1)
        layout.addWidget(opt_strict)
        layout.addWidget(opt_loose)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        if dlg.exec() != QDialog.Accepted:
            return None
        return 1.0e-3 if opt_strict.isChecked() else 5.0e-2

    def _prompt_symmetrize_coords(self) -> bool | None:
        dlg = QDialog(self)
        dlg.setWindowTitle("Symmetrize coordinates")
        layout = QVBoxLayout(dlg)
        title = QLabel("Symmetrize coordinates (replace xyzin geometry)?")
        title.setStyleSheet("font-weight: 600;")
        layout.addWidget(title)

        group = QButtonGroup(dlg)
        group.setExclusive(True)
        opt_yes = QRadioButton("Yes (apply symmetrized coordinates to xyzin)", dlg)
        opt_no = QRadioButton("No (report only)", dlg)
        opt_no.setChecked(True)
        group.addButton(opt_yes, 0)
        group.addButton(opt_no, 1)
        layout.addWidget(opt_yes)
        layout.addWidget(opt_no)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        if dlg.exec() != QDialog.Accepted:
            return None
        return bool(opt_yes.isChecked())

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
            self.basic_representation,
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
            self.basic_representation,
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

    def _open_isotopologues_dialog(self):
        if not self.xyzin_path.exists():
            QMessageBox.warning(self, "Isotopologues", "No xyzin available.")
            return
        try:
            _nat, _comment, symbols, _coords, _tail = read_xyz_full_from_xyzin(
                str(self.xyzin_path)
            )
        except Exception:
            QMessageBox.warning(self, "Isotopologues", "Failed to read xyzin.")
            return

        default_isos = []
        for s in symbols:
            z = atomic_number(s)
            iso = get_default_isotope(int(z))
            default_isos.append(iso.A if iso is not None else None)

        dlg = QDialog(self)
        dlg.setWindowTitle("Isotopologues")
        dlg.resize(620, 420)
        layout = QVBoxLayout(dlg)

        title = QLabel("Define isotopologues (only non-default isotopes).")
        title.setStyleSheet("font-weight: 600;")
        layout.addWidget(title)

        row = QHBoxLayout()
        row.addWidget(QLabel("Number of isotopologues:"))
        iso_count = QSpinBox(dlg)
        iso_count.setRange(1, 100)
        iso_count.setValue(1)
        row.addWidget(iso_count)
        row.addStretch()
        layout.addLayout(row)

        table = QTableWidget(dlg)
        table.setColumnCount(3)
        table.setHorizontalHeaderLabels(["Isotopologue #", "Atom (1-based)", "Isotope mass A"])
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        table.setRowCount(0)
        layout.addWidget(table)

        btn_row = QHBoxLayout()
        add_row = QPushButton("Add row")
        del_row = QPushButton("Remove row")
        btn_row.addWidget(add_row)
        btn_row.addWidget(del_row)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        err_label = QLabel("")
        err_label.setStyleSheet("color: #b00020;")
        layout.addWidget(err_label)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        def _build_atom_combo():
            combo = QComboBox(dlg)
            for i, s in enumerate(symbols, start=1):
                d = default_isos[i - 1]
                combo.addItem(f"{i}: {s} (A={d})", i)
            return combo

        def _build_iso_combo(atom_idx: int):
            combo = QComboBox(dlg)
            s = symbols[atom_idx - 1]
            z = atomic_number(s)
            default_a = default_isos[atom_idx - 1]
            for iso in get_isotopes(int(z)) or []:
                if default_a is not None and int(iso.A) == int(default_a):
                    continue
                combo.addItem(f"{iso.A}", int(iso.A))
            return combo

        def _set_row_widgets(row: int):
            iso_spin = QSpinBox(dlg)
            iso_spin.setRange(1, 100)
            iso_spin.setValue(1)
            atom_combo = _build_atom_combo()
            iso_combo = _build_iso_combo(1)

            def _atom_changed():
                idx = int(atom_combo.currentData() or 1)
                iso_combo.clear()
                s = symbols[idx - 1]
                z = atomic_number(s)
                default_a = default_isos[idx - 1]
                for iso in get_isotopes(int(z)) or []:
                    if default_a is not None and int(iso.A) == int(default_a):
                        continue
                    iso_combo.addItem(f"{iso.A}", int(iso.A))

            atom_combo.currentIndexChanged.connect(_atom_changed)
            table.setCellWidget(row, 0, iso_spin)
            table.setCellWidget(row, 1, atom_combo)
            table.setCellWidget(row, 2, iso_combo)

        def _add_row():
            row = table.rowCount()
            table.setRowCount(row + 1)
            _set_row_widgets(row)

        def _del_row():
            if table.rowCount() > 0:
                table.setRowCount(table.rowCount() - 1)

        add_row.clicked.connect(_add_row)
        del_row.clicked.connect(_del_row)
        _add_row()

        if dlg.exec() != QDialog.Accepted:
            return

        n_iso = int(iso_count.value())
        rows = []
        for r in range(table.rowCount()):
            w0 = table.cellWidget(r, 0)
            w1 = table.cellWidget(r, 1)
            w2 = table.cellWidget(r, 2)
            if w0 is None or w1 is None or w2 is None:
                continue
            i_iso = int(w0.value())
            atom_idx = int(w1.currentData() or 0)
            mass_a = int(w2.currentData() or 0)
            rows.append((i_iso, atom_idx, mass_a))

        if not rows:
            err_label.setText("No valid substitutions specified.")
            return

        # Validate ranges and defaults
        for i_iso, atom_idx, mass_a in rows:
            if i_iso < 1 or i_iso > n_iso:
                err_label.setText(f"Isotopologue # out of range: {i_iso}")
                return
            if atom_idx < 1 or atom_idx > len(symbols):
                err_label.setText(f"Atom index out of range: {atom_idx}")
                return
            default_a = default_isos[atom_idx - 1]
            if default_a is not None and int(mass_a) == int(default_a):
                err_label.setText(
                    f"Atom {atom_idx} uses default isotope {default_a}; choose a non-default A."
                )
                return
            # Validate isotope exists for element
            z = atomic_number(symbols[atom_idx - 1])
            allowed = {int(iso.A) for iso in (get_isotopes(int(z)) or [])}
            if int(mass_a) not in allowed:
                err_label.setText(
                    f"Isotope A={mass_a} not available for element {symbols[atom_idx - 1]}."
                )
                return

        # Build output
        by_iso = {i: [] for i in range(1, n_iso + 1)}
        for i_iso, atom_idx, mass_a in rows:
            by_iso[i_iso].append((atom_idx, mass_a))

        out_lines = [str(n_iso)]
        for i in range(1, n_iso + 1):
            subs = by_iso.get(i, [])
            out_lines.append(str(len(subs)))
            for atom_idx, mass_a in subs:
                out_lines.append(f"{atom_idx} {mass_a}")

        out_path = self.working_dir / "isotopologues.txt"
        out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
        records = []
        for i in range(1, n_iso + 1):
            records.append(
                XyzinIsotopologueRecord(
                    label=f"iso_{i:03d}",
                    substitutions={atom_idx: mass_a for atom_idx, mass_a in by_iso.get(i, [])},
                )
            )
        merge_xyzin_isotopologue_records(self.xyzin_path, tuple(records))
        QMessageBox.information(
            self,
            "Isotopologues",
            f"Saved: {out_path}\nUpdated: {self.xyzin_path} (#ISOTOPOLOGUES)",
        )

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
