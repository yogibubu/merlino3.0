from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from merlino_vpt2_vci import (
    VCIOptions,
    load_force_field,
    run_gf_report_from_fchk,
    run_vpt2_vci_report,
)


class VPT2VCIWindow(QMainWindow):
    """Dedicated GUI for GF/PED and anharmonic VPT2/VCI workflows."""

    def __init__(self, workdir: Path, repo_root: Path | None = None, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.repo_root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]

        self.setWindowTitle("Merlino GF / VPT2-VCI")
        self.resize(980, 780)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        header = QLabel("GF / VPT2-VCI")
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        note = QLabel(
            "Gaussian/FCHK is only an input adapter here. GF, PED, VPT2 and VCI "
            "run on canonical Merlino data structures and Merlino non-redundant GICs."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        input_group = QGroupBox("Inputs")
        input_layout = QVBoxLayout(input_group)

        row_fchk = QHBoxLayout()
        row_fchk.addWidget(QLabel("FCHK Hessian/QFF:"))
        self.fchk_edit = QLineEdit(str(self.workdir / "gauin.fchk"))
        row_fchk.addWidget(self.fchk_edit)
        browse_fchk = QPushButton("Browse")
        browse_fchk.clicked.connect(self._browse_fchk)
        row_fchk.addWidget(browse_fchk)
        input_layout.addLayout(row_fchk)

        row_qff = QHBoxLayout()
        row_qff.addWidget(QLabel("Indexed QFF text:"))
        self.qff_edit = QLineEdit("")
        row_qff.addWidget(self.qff_edit)
        browse_qff = QPushButton("Browse")
        browse_qff.clicked.connect(self._browse_qff)
        row_qff.addWidget(browse_qff)
        input_layout.addLayout(row_qff)

        row_latest = QHBoxLayout()
        latest_fchk = QPushButton("Use Latest FCHK")
        latest_fchk.clicked.connect(lambda: self._select_latest_fchk())
        row_latest.addWidget(latest_fchk)
        row_latest.addStretch()
        input_layout.addLayout(row_latest)
        layout.addWidget(input_group)

        gf_group = QGroupBox("GF / PED")
        gf_layout = QVBoxLayout(gf_group)
        gf_note = QLabel(
            "Reads Cartesian Hessian and geometry from FCHK, builds Merlino GICs, "
            "then solves Wilson GF and reports PED in non-redundant GICs."
        )
        gf_note.setWordWrap(True)
        gf_layout.addWidget(gf_note)
        self.run_gf_button = QPushButton("Run GF / PED")
        self.run_gf_button.clicked.connect(self.run_gf)
        gf_layout.addWidget(self.run_gf_button)
        layout.addWidget(gf_group)

        vci_group = QGroupBox("VPT2 / VCI")
        vci_layout = QVBoxLayout(vci_group)

        row_basis = QHBoxLayout()
        row_basis.addWidget(QLabel("Max total quanta:"))
        self.max_quanta_edit = QLineEdit("2")
        row_basis.addWidget(self.max_quanta_edit)
        row_basis.addWidget(QLabel("Roots:"))
        self.roots_edit = QLineEdit("6")
        row_basis.addWidget(self.roots_edit)
        row_basis.addWidget(QLabel("Active modes (1-based):"))
        self.active_modes_edit = QLineEdit("")
        self.active_modes_edit.setPlaceholderText("blank = all, e.g. 1,2,5")
        row_basis.addWidget(self.active_modes_edit)
        vci_layout.addLayout(row_basis)

        row_prune = QHBoxLayout()
        row_prune.addWidget(QLabel("Freq min cm-1:"))
        self.freq_min_edit = QLineEdit("")
        row_prune.addWidget(self.freq_min_edit)
        row_prune.addWidget(QLabel("Freq max cm-1:"))
        self.freq_max_edit = QLineEdit("")
        row_prune.addWidget(self.freq_max_edit)
        row_prune.addWidget(QLabel("Basis cutoff cm-1:"))
        self.basis_cutoff_edit = QLineEdit("")
        row_prune.addWidget(self.basis_cutoff_edit)
        row_prune.addWidget(QLabel("Force threshold cm-1:"))
        self.force_threshold_edit = QLineEdit("0.0")
        row_prune.addWidget(self.force_threshold_edit)
        vci_layout.addLayout(row_prune)

        row_limits = QHBoxLayout()
        row_limits.addWidget(QLabel("Mode max quanta:"))
        self.mode_max_edit = QLineEdit("")
        self.mode_max_edit.setPlaceholderText("blank or per original mode, e.g. 3,3,2")
        row_limits.addWidget(self.mode_max_edit)
        row_limits.addWidget(QLabel("Class limits 1/2/3/4-mode:"))
        self.class_limits_edit = QLineEdit("")
        self.class_limits_edit.setPlaceholderText("e.g. 1:1-2;2:2-3;3:3-4;4:4-4")
        row_limits.addWidget(self.class_limits_edit)
        vci_layout.addLayout(row_limits)

        row_run = QHBoxLayout()
        self.run_vpt2_vci_button = QPushButton("Run VPT2 / VCI")
        self.run_vpt2_vci_button.clicked.connect(self.run_vpt2_vci)
        row_run.addWidget(self.run_vpt2_vci_button)
        row_run.addStretch()
        vci_layout.addLayout(row_run)
        layout.addWidget(vci_group)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setPlaceholderText("GF/PED and VPT2/VCI reports will appear here.")

        actions = QHBoxLayout()
        clear_button = QPushButton("Clear Output")
        clear_button.clicked.connect(self.output_text.clear)
        actions.addWidget(clear_button)
        export_button = QPushButton("Export Report")
        export_button.clicked.connect(self.export_report)
        actions.addWidget(export_button)
        save_preset_button = QPushButton("Save Preset")
        save_preset_button.clicked.connect(self.save_preset)
        actions.addWidget(save_preset_button)
        load_preset_button = QPushButton("Load Preset")
        load_preset_button.clicked.connect(self.load_preset)
        actions.addWidget(load_preset_button)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        actions.addWidget(close_button)
        actions.addStretch()
        layout.addLayout(actions)

        layout.addWidget(self.output_text)

    def run_gf(self, *, show_message: bool = True) -> None:
        try:
            fchk_path = self._required_existing_path(self.fchk_edit.text(), "FCHK")
            report = run_gf_report_from_fchk(fchk_path)
            self.output_text.setPlainText(report.text)
        except Exception as exc:
            self._fail("GF / PED failed", exc, show_message)

    def run_vpt2_vci(self, *, show_message: bool = True) -> None:
        try:
            qff = self._load_force_field()
            max_quanta = self._parse_int(self.max_quanta_edit.text(), "Max total quanta", minimum=0)
            roots = self._parse_int(self.roots_edit.text(), "Roots", minimum=1)
            report = run_vpt2_vci_report(qff, max_quanta=max_quanta, roots=roots, options=self._vci_options())
            self.output_text.setPlainText(report.text)
        except Exception as exc:
            self._fail("VPT2 / VCI failed", exc, show_message)

    def _load_force_field(self):
        qff_path = self._optional_existing_path(self.qff_edit.text())
        fchk_path = self._optional_existing_path(self.fchk_edit.text())
        return load_force_field(fchk_path=fchk_path, qff_path=qff_path)

    def _vci_options(self) -> VCIOptions:
        return VCIOptions(
            active_modes=self._parse_active_modes(),
            frequency_min_cm=self._parse_optional_float(self.freq_min_edit.text(), "Freq min cm-1"),
            frequency_max_cm=self._parse_optional_float(self.freq_max_edit.text(), "Freq max cm-1"),
            basis_energy_cutoff_cm=self._parse_optional_float(self.basis_cutoff_edit.text(), "Basis cutoff cm-1"),
            mode_max_quanta=self._parse_optional_int_tuple(self.mode_max_edit.text(), "Mode max quanta"),
            excitation_class_limits=self._parse_class_limits(),
            force_constant_threshold_cm=self._parse_optional_float(
                self.force_threshold_edit.text(), "Force threshold cm-1", default=0.0
            ),
        )

    def export_report(self, path: Path | None = None, *, show_message: bool = True) -> Path | None:
        text = self.output_text.toPlainText()
        if not text.strip():
            if show_message:
                QMessageBox.warning(self, "GF / VPT2-VCI", "No report to export.")
            return None
        if path is None:
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Export GF / VPT2-VCI report",
                str(self.workdir / "vpt2_vci_report.txt"),
                "Text files (*.txt);;All files (*)",
            )
            if not selected:
                return None
            path = Path(selected)
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
        if show_message:
            QMessageBox.information(self, "GF / VPT2-VCI", f"Report written: {path}")
        return path

    def save_preset(self, path: Path | None = None, *, show_message: bool = True) -> Path | None:
        if path is None:
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Save VPT2/VCI preset",
                str(self.workdir / "vpt2_vci_preset.json"),
                "JSON files (*.json);;All files (*)",
            )
            if not selected:
                return None
            path = Path(selected)
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self._settings_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if show_message:
            QMessageBox.information(self, "GF / VPT2-VCI", f"Preset written: {path}")
        return path

    def load_preset(self, path: Path | None = None, *, show_message: bool = True) -> Path | None:
        if path is None:
            selected, _ = QFileDialog.getOpenFileName(
                self,
                "Load VPT2/VCI preset",
                str(self.workdir),
                "JSON files (*.json);;All files (*)",
            )
            if not selected:
                return None
            path = Path(selected)
        path = Path(path).expanduser()
        data = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("Preset root must be a JSON object")
        self._apply_settings_dict(data)
        if show_message:
            QMessageBox.information(self, "GF / VPT2-VCI", f"Preset loaded: {path}")
        return path

    def _settings_dict(self) -> dict[str, str]:
        return {
            "fchk_path": self.fchk_edit.text(),
            "qff_path": self.qff_edit.text(),
            "max_quanta": self.max_quanta_edit.text(),
            "roots": self.roots_edit.text(),
            "active_modes": self.active_modes_edit.text(),
            "frequency_min_cm": self.freq_min_edit.text(),
            "frequency_max_cm": self.freq_max_edit.text(),
            "basis_cutoff_cm": self.basis_cutoff_edit.text(),
            "force_threshold_cm": self.force_threshold_edit.text(),
            "mode_max_quanta": self.mode_max_edit.text(),
            "class_limits": self.class_limits_edit.text(),
        }

    def _apply_settings_dict(self, data: dict[str, object]) -> None:
        mapping = {
            "fchk_path": self.fchk_edit,
            "qff_path": self.qff_edit,
            "max_quanta": self.max_quanta_edit,
            "roots": self.roots_edit,
            "active_modes": self.active_modes_edit,
            "frequency_min_cm": self.freq_min_edit,
            "frequency_max_cm": self.freq_max_edit,
            "basis_cutoff_cm": self.basis_cutoff_edit,
            "force_threshold_cm": self.force_threshold_edit,
            "mode_max_quanta": self.mode_max_edit,
            "class_limits": self.class_limits_edit,
        }
        for key, widget in mapping.items():
            if key in data:
                widget.setText("" if data[key] is None else str(data[key]))

    def _browse_fchk(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select FCHK", str(self.workdir), "FCHK files (*.fchk *.fch);;All files (*)")
        if path:
            self.fchk_edit.setText(path)

    def _browse_qff(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select indexed QFF", str(self.workdir), "QFF text (*.qff *.txt);;All files (*)")
        if path:
            self.qff_edit.setText(path)

    def _select_latest_fchk(self, *, show_message: bool = True) -> None:
        candidates = sorted(
            list(self.workdir.glob("*.fchk")) + list(self.workdir.glob("*.fch")),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if not candidates:
            if show_message:
                QMessageBox.warning(self, "FCHK", f"No FCHK files found in {self.workdir}")
            return
        self.fchk_edit.setText(str(candidates[0]))

    def _required_existing_path(self, raw: str, label: str) -> Path:
        path = Path(raw.strip()).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
        return path

    def _optional_existing_path(self, raw: str) -> Path | None:
        text = raw.strip()
        if not text:
            return None
        path = Path(text).expanduser()
        if not path.exists():
            return None
        return path

    def _parse_active_modes(self) -> tuple[int, ...] | None:
        values = self._parse_optional_int_tuple(self.active_modes_edit.text(), "Active modes")
        if values is None:
            return None
        if any(value < 1 for value in values):
            raise ValueError("Active modes are one-based and must be positive")
        return tuple(value - 1 for value in values)

    def _parse_class_limits(self) -> dict[int, tuple[int, int | None]]:
        text = self.class_limits_edit.text().strip()
        if not text:
            return {}
        limits: dict[int, tuple[int, int | None]] = {}
        for item in text.replace(",", ";").split(";"):
            item = item.strip()
            if not item:
                continue
            if ":" not in item:
                raise ValueError("Class limits must use n:min-max records")
            key_raw, value_raw = item.split(":", 1)
            n_modes = self._parse_int(key_raw, "Class limit key", minimum=1)
            if "-" in value_raw:
                lo_raw, hi_raw = value_raw.split("-", 1)
                qmin = self._parse_int(lo_raw, "Class limit minimum", minimum=0)
                qmax = None if hi_raw.strip() == "*" else self._parse_int(hi_raw, "Class limit maximum", minimum=qmin)
            else:
                qmin = self._parse_int(value_raw, "Class limit minimum", minimum=0)
                qmax = qmin
            limits[n_modes] = (qmin, qmax)
        return limits

    @staticmethod
    def _parse_optional_int_tuple(raw: str, label: str) -> tuple[int, ...] | None:
        text = raw.strip()
        if not text:
            return None
        values = []
        for token in text.replace(";", ",").split(","):
            token = token.strip()
            if token:
                values.append(VPT2VCIWindow._parse_int(token, label, minimum=0))
        return tuple(values)

    @staticmethod
    def _parse_optional_float(raw: str, label: str, *, default: float | None = None) -> float | None:
        text = raw.strip()
        if not text:
            return default
        try:
            return float(text)
        except ValueError as exc:
            raise ValueError(f"{label} must be numeric") from exc

    @staticmethod
    def _parse_int(raw: str, label: str, *, minimum: int) -> int:
        try:
            value = int(raw.strip())
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer") from exc
        if value < minimum:
            raise ValueError(f"{label} must be >= {minimum}")
        return value

    def _fail(self, title: str, exc: Exception, show_message: bool) -> None:
        message = f"{title}: {exc}"
        self.output_text.setPlainText(message)
        if show_message:
            QMessageBox.critical(self, title, str(exc))
