from __future__ import annotations

import json
from pathlib import Path

from merlino_core import build_run_manifest
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

from merlino_gf import (
    run_gic_gf_report_from_fchk,
    run_gf_report_from_fchk,
    write_csv_tables,
)


class GFWindow(QMainWindow):
    """Dedicated GUI for Wilson GF/PED workflows."""

    def __init__(self, workdir: Path, repo_root: Path | None = None, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.repo_root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]
        self.last_report = None

        self.setWindowTitle("Merlino GF / PED")
        self.resize(980, 700)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        header = QLabel("GF / PED in Merlino GICs")
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        note = QLabel(
            "This branch reads a Cartesian Hessian, evaluates Merlino GICs and "
            "B matrices, transforms the Hessian to internal coordinates, and "
            "solves Wilson GF. It is separate from VPT2/VCI, which works in "
            "Cartesian normal modes."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        input_group = QGroupBox("Inputs")
        input_layout = QVBoxLayout(input_group)

        row_fchk = QHBoxLayout()
        row_fchk.addWidget(QLabel("FCHK Hessian:"))
        self.fchk_edit = QLineEdit(str(self.workdir / "gauin.fchk"))
        row_fchk.addWidget(self.fchk_edit)
        browse_fchk = QPushButton("Browse")
        browse_fchk.clicked.connect(self._browse_fchk)
        row_fchk.addWidget(browse_fchk)
        input_layout.addLayout(row_fchk)

        row_gic = QHBoxLayout()
        row_gic.addWidget(QLabel("Frozen GIC definition:"))
        self.gic_schema_edit = QLineEdit("")
        self.gic_schema_edit.setPlaceholderText("optional: gic_definition.json from gic-define")
        row_gic.addWidget(self.gic_schema_edit)
        browse_gic = QPushButton("Browse")
        browse_gic.clicked.connect(self._browse_gic_schema)
        row_gic.addWidget(browse_gic)
        input_layout.addLayout(row_gic)

        row_gic_geom = QHBoxLayout()
        row_gic_geom.addWidget(QLabel("B geometry override:"))
        self.gic_geometry_edit = QLineEdit("")
        self.gic_geometry_edit.setPlaceholderText("optional XYZ/COM/GJF; blank uses FCHK geometry")
        row_gic_geom.addWidget(self.gic_geometry_edit)
        browse_gic_geom = QPushButton("Browse")
        browse_gic_geom.clicked.connect(self._browse_gic_geometry)
        row_gic_geom.addWidget(browse_gic_geom)
        input_layout.addLayout(row_gic_geom)

        row_scale = QHBoxLayout()
        row_scale.addWidget(QLabel("Pulay scale file:"))
        self.gic_scale_edit = QLineEdit("")
        self.gic_scale_edit.setPlaceholderText("optional: GIC001 0.98; default 0.97")
        row_scale.addWidget(self.gic_scale_edit)
        browse_scale = QPushButton("Browse")
        browse_scale.clicked.connect(self._browse_gic_scale)
        row_scale.addWidget(browse_scale)
        input_layout.addLayout(row_scale)

        row_latest = QHBoxLayout()
        latest_fchk = QPushButton("Use Latest FCHK")
        latest_fchk.clicked.connect(lambda: self._select_latest_fchk())
        row_latest.addWidget(latest_fchk)
        row_latest.addStretch()
        input_layout.addLayout(row_latest)
        layout.addWidget(input_group)

        run_group = QGroupBox("GF / PED")
        run_layout = QVBoxLayout(run_group)
        self.run_gf_button = QPushButton("Run GF / PED")
        self.run_gf_button.clicked.connect(self.run_gf)
        run_layout.addWidget(self.run_gf_button)
        layout.addWidget(run_group)

        actions = QHBoxLayout()
        clear_button = QPushButton("Clear Output")
        clear_button.clicked.connect(lambda: self.output_text_clear())
        actions.addWidget(clear_button)
        export_button = QPushButton("Export Report")
        export_button.clicked.connect(lambda: self.export_report())
        actions.addWidget(export_button)
        export_csv_button = QPushButton("Export CSVs")
        export_csv_button.clicked.connect(lambda: self.export_csvs())
        actions.addWidget(export_csv_button)
        save_preset_button = QPushButton("Save Preset")
        save_preset_button.clicked.connect(lambda: self.save_preset())
        actions.addWidget(save_preset_button)
        load_preset_button = QPushButton("Load Preset")
        load_preset_button.clicked.connect(lambda: self.load_preset())
        actions.addWidget(load_preset_button)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        actions.addWidget(close_button)
        actions.addStretch()
        layout.addLayout(actions)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setPlaceholderText("GF/PED report will appear here.")
        layout.addWidget(self.output_text)

    def output_text_clear(self) -> None:
        self.output_text.clear()

    def run_gf(self, *, show_message: bool = True) -> None:
        try:
            fchk_path = self._required_existing_path(self.fchk_edit.text(), "FCHK")
            gic_schema = self._optional_existing_path(self.gic_schema_edit.text())
            if gic_schema is None:
                report = run_gf_report_from_fchk(fchk_path)
            else:
                report = run_gic_gf_report_from_fchk(
                    fchk_path,
                    gic_schema,
                    geometry_path=self._optional_existing_path(self.gic_geometry_edit.text()),
                    scale_path=self._optional_existing_path(self.gic_scale_edit.text()),
                )
            self.last_report = report
            self.output_text.setPlainText(report.text)
            self._write_gf_manifest(fchk_path)
        except Exception as exc:
            self._fail("GF / PED failed", exc, show_message)

    def export_report(self, path: Path | None = None, *, show_message: bool = True) -> Path | None:
        text = self.output_text.toPlainText()
        if not text.strip():
            if show_message:
                QMessageBox.warning(self, "GF / PED", "No report to export.")
            return None
        if path is None:
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Export GF/PED report",
                str(self.workdir / "gf_ped_report.txt"),
                "Text files (*.txt);;All files (*)",
            )
            if not selected:
                return None
            path = Path(selected)
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
        if show_message:
            QMessageBox.information(self, "GF / PED", f"Report written: {path}")
        return path

    def export_csvs(self, outdir: Path | None = None, *, show_message: bool = True) -> dict[str, Path]:
        if self.last_report is None:
            if show_message:
                QMessageBox.warning(self, "GF / PED", "No completed GF/PED result to export.")
            return {}
        if outdir is None:
            selected = QFileDialog.getExistingDirectory(self, "Export GF/PED CSV tables", str(self.workdir))
            if not selected:
                return {}
            outdir = Path(selected)
        written = write_csv_tables(self.last_report, Path(outdir))
        if show_message:
            QMessageBox.information(self, "GF / PED", f"CSV tables written: {Path(outdir)}")
        return written

    def _write_gf_manifest(self, fchk_path: Path) -> Path:
        inputs = {"fchk": fchk_path}
        gic_schema = self._optional_existing_path(self.gic_schema_edit.text())
        gic_geometry = self._optional_existing_path(self.gic_geometry_edit.text())
        scale_path = self._optional_existing_path(self.gic_scale_edit.text())
        if gic_schema is not None:
            inputs["gic_definition"] = gic_schema
        if gic_geometry is not None:
            inputs["geometry"] = gic_geometry
        if scale_path is not None:
            inputs["scale_file"] = scale_path
        workflow = "gic_gf" if gic_schema is not None else "gf"
        return build_run_manifest(
            workflow=workflow,
            status="completed",
            run_dir=self.workdir,
            inputs=inputs,
            backend={
                "adapter": "gaussian-fchk",
                "solver": "python",
                "gui": "advanced.gf_window",
                "coordinate_model": "frozen-gic-definition" if gic_schema is not None else "generated-merlino-gics",
            },
        ).write(self.workdir / ("gic_gf_manifest.json" if workflow == "gic_gf" else "gf_manifest.json"))

    def save_preset(self, path: Path | None = None, *, show_message: bool = True) -> Path | None:
        if path is None:
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Save GF/PED preset",
                str(self.workdir / "gf_ped_preset.json"),
                "JSON files (*.json);;All files (*)",
            )
            if not selected:
                return None
            path = Path(selected)
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self._settings_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if show_message:
            QMessageBox.information(self, "GF / PED", f"Preset written: {path}")
        return path

    def load_preset(self, path: Path | None = None, *, show_message: bool = True) -> Path | None:
        if path is None:
            selected, _ = QFileDialog.getOpenFileName(
                self,
                "Load GF/PED preset",
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
            QMessageBox.information(self, "GF / PED", f"Preset loaded: {path}")
        return path

    def _settings_dict(self) -> dict[str, str]:
        return {
            "fchk_path": self.fchk_edit.text(),
            "gic_schema_path": self.gic_schema_edit.text(),
            "gic_geometry_path": self.gic_geometry_edit.text(),
            "gic_scale_path": self.gic_scale_edit.text(),
        }

    def _apply_settings_dict(self, data: dict[str, object]) -> None:
        mapping = {
            "fchk_path": self.fchk_edit,
            "gic_schema_path": self.gic_schema_edit,
            "gic_geometry_path": self.gic_geometry_edit,
            "gic_scale_path": self.gic_scale_edit,
        }
        for key, widget in mapping.items():
            if key in data:
                widget.setText("" if data[key] is None else str(data[key]))

    def _browse_fchk(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select FCHK", str(self.workdir), "FCHK files (*.fchk *.fch);;All files (*)")
        if path:
            self.fchk_edit.setText(path)

    def _browse_gic_schema(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select GIC definition", str(self.workdir), "JSON files (*.json);;All files (*)")
        if path:
            self.gic_schema_edit.setText(path)

    def _browse_gic_geometry(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select B-matrix geometry",
            str(self.workdir),
            "Cartesian geometry (*.xyz *.com *.gjf);;All files (*)",
        )
        if path:
            self.gic_geometry_edit.setText(path)

    def _browse_gic_scale(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Pulay scaling file",
            str(self.workdir),
            "Scaling files (*.csv *.txt *.scale);;All files (*)",
        )
        if path:
            self.gic_scale_edit.setText(path)

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

    def _fail(self, title: str, exc: Exception, show_message: bool) -> None:
        message = f"{title}: {exc}"
        self.output_text.setPlainText(message)
        if show_message:
            QMessageBox.critical(self, title, str(exc))
