from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QDoubleSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from merlino_fit.survibfit.fragment_pipeline import run_fragment_pipeline


class FragmentPipelineWindow(QDialog):
    def __init__(self, working_dir: Path, parent=None):
        super().__init__(parent)
        self.working_dir = Path(working_dir)
        self._last_report = None

        self.setWindowTitle("Fragment Pipeline (MVP)")
        self.resize(860, 620)

        root = QVBoxLayout(self)
        title = QLabel(
            "Fragment similarity + GAP detection (SE preferred over PCS2 on duplicate molecules)"
        )
        title.setStyleSheet("font-weight: 600;")
        root.addWidget(title)

        panel = QWidget(self)
        form = QFormLayout(panel)
        root.addWidget(panel)

        self.query_xyz = QLineEdit(panel)
        self.se_dir = QLineEdit(panel)
        self.pcs2_dir = QLineEdit(panel)
        self.out_dir = QLineEdit(panel)
        self.top_k = QSpinBox(panel)
        self.top_k.setRange(1, 200)
        self.top_k.setValue(5)
        self.gap_threshold = QDoubleSpinBox(panel)
        self.gap_threshold.setRange(0.0, 1.0)
        self.gap_threshold.setDecimals(3)
        self.gap_threshold.setValue(0.75)

        row_q = QHBoxLayout()
        row_q.addWidget(self.query_xyz)
        btn_q = QPushButton("Browse…", panel)
        row_q.addWidget(btn_q)
        form.addRow("Query XYZ", row_q)

        row_se = QHBoxLayout()
        row_se.addWidget(self.se_dir)
        btn_se = QPushButton("Browse…", panel)
        row_se.addWidget(btn_se)
        form.addRow("SE library dir", row_se)

        row_pcs2 = QHBoxLayout()
        row_pcs2.addWidget(self.pcs2_dir)
        btn_pcs2 = QPushButton("Browse…", panel)
        row_pcs2.addWidget(btn_pcs2)
        form.addRow("PCS2 library dir", row_pcs2)

        row_out = QHBoxLayout()
        row_out.addWidget(self.out_dir)
        btn_out = QPushButton("Browse…", panel)
        row_out.addWidget(btn_out)
        form.addRow("Output dir", row_out)

        form.addRow("Top-k per fragment", self.top_k)
        form.addRow("GAP threshold", self.gap_threshold)

        run_btn = QPushButton("Run fragment pipeline", panel)
        form.addRow("", run_btn)

        self.output = QTextEdit(self)
        self.output.setReadOnly(True)
        root.addWidget(self.output)

        footer = QHBoxLayout()
        root.addLayout(footer)
        save_btn = QPushButton("Save JSON report")
        close_btn = QPushButton("Close")
        footer.addWidget(save_btn)
        footer.addStretch()
        footer.addWidget(close_btn)

        btn_q.clicked.connect(lambda: self._pick_file(self.query_xyz))
        btn_se.clicked.connect(lambda: self._pick_dir(self.se_dir))
        btn_pcs2.clicked.connect(lambda: self._pick_dir(self.pcs2_dir))
        btn_out.clicked.connect(lambda: self._pick_dir(self.out_dir))
        run_btn.clicked.connect(self._run_pipeline)
        save_btn.clicked.connect(self._save_json_report)
        close_btn.clicked.connect(self.accept)

        # Defaults
        self.se_dir.setText(str(self.working_dir.parent / "projects" / "se_library"))
        self.pcs2_dir.setText(str(self.working_dir.parent / "projects" / "pcs2_library"))
        self.out_dir.setText(str(self.working_dir / "fragment_reports"))

    def _pick_file(self, target: QLineEdit):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select XYZ file", str(self.working_dir), "XYZ files (*.xyz);;All files (*)"
        )
        if path:
            target.setText(path)

    def _pick_dir(self, target: QLineEdit):
        path = QFileDialog.getExistingDirectory(self, "Select directory", str(self.working_dir))
        if path:
            target.setText(path)

    def _run_pipeline(self):
        query = Path(self.query_xyz.text().strip())
        se_dir = Path(self.se_dir.text().strip())
        pcs2_dir = Path(self.pcs2_dir.text().strip())
        out_dir = Path(self.out_dir.text().strip())
        if not query.exists():
            QMessageBox.warning(self, "Fragment pipeline", "Please select a valid query XYZ file.")
            return
        if not se_dir.is_dir() or not pcs2_dir.is_dir():
            QMessageBox.warning(
                self, "Fragment pipeline", "Please select valid SE/PCS2 library directories."
            )
            return
        try:
            report = run_fragment_pipeline(
                query,
                se_dir,
                pcs2_dir,
                out_dir,
                top_k=int(self.top_k.value()),
                gap_threshold=float(self.gap_threshold.value()),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Fragment pipeline error", str(exc))
            return

        self._last_report = report
        lines = [
            "Fragment pipeline completed",
            f"Query: {report['query_xyz']}",
            f"Query fragments: {report['query_fragments']}",
            f"Library molecules: {report['library_xyz_count']}",
            f"Library fragments: {report['library_fragments']}",
            f"GAP threshold: {report['gap_threshold']:.3f}",
            f"GAP fragments: {len(report['to_curate'])}",
            "",
            "Top fragment statuses:",
        ]
        for row in report["fragments"][:10]:
            lines.append(
                f" - {row['query_fragment_id']}: best={row['best_similarity']:.6f} "
                f"status={row['status']}"
            )
        lines += [
            "",
            f"Saved: {out_dir / 'fragment_pipeline.json'}",
            f"Saved: {out_dir / 'to_curate.json'}",
            f"Saved: {out_dir / 'fragment_pipeline.md'}",
        ]
        self.output.setPlainText("\n".join(lines))

    def _save_json_report(self):
        if self._last_report is None:
            QMessageBox.information(self, "Fragment pipeline", "No result available yet.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save fragment pipeline JSON",
            str(self.working_dir / "fragment_pipeline.json"),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            Path(path).write_text(json.dumps(self._last_report, indent=2) + "\n", encoding="utf-8")
        except Exception as exc:
            QMessageBox.critical(self, "Save error", str(exc))
            return
        QMessageBox.information(self, "Fragment pipeline", f"Saved report:\n{path}")
