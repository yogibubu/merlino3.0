from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
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
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from merlino_fit.survibfit.synthon_similarity import (
    compare_against_library,
    compare_molecules,
)


class SimilarityWindow(QDialog):
    def __init__(self, working_dir: Path, parent=None):
        super().__init__(parent)
        self.working_dir = Path(working_dir)
        self._last_result = None

        self.setWindowTitle("Synthon Similarity")
        self.resize(860, 620)

        root = QVBoxLayout(self)
        title = QLabel(
            "Molecule similarity using synthon Gaussian models "
            "(charge, covalency, delocalization, strain, Zeff)"
        )
        title.setStyleSheet("font-weight: 600;")
        root.addWidget(title)

        self.tabs = QTabWidget(self)
        root.addWidget(self.tabs)

        self._build_pair_tab()
        self._build_library_tab()

        self.output = QTextEdit(self)
        self.output.setReadOnly(True)
        root.addWidget(self.output)

        footer = QHBoxLayout()
        root.addLayout(footer)
        self.save_json_btn = QPushButton("Save JSON report")
        self.close_btn = QPushButton("Close")
        footer.addWidget(self.save_json_btn)
        footer.addStretch()
        footer.addWidget(self.close_btn)

        self.save_json_btn.clicked.connect(self._save_json_report)
        self.close_btn.clicked.connect(self.accept)

    def _build_pair_tab(self):
        tab = QWidget(self)
        form = QFormLayout(tab)

        self.pair_a = QLineEdit(tab)
        self.pair_b = QLineEdit(tab)
        self.pair_cov = QComboBox(tab)
        self.pair_cov.addItems(["full", "diag"])
        self.pair_reg = QDoubleSpinBox(tab)
        self.pair_reg.setRange(1.0e-8, 1.0)
        self.pair_reg.setDecimals(6)
        self.pair_reg.setValue(5.0e-2)
        self.pair_std = QCheckBox("Standardize features", tab)
        self.pair_std.setChecked(True)

        row_a = QHBoxLayout()
        row_a.addWidget(self.pair_a)
        btn_a = QPushButton("Browse…", tab)
        row_a.addWidget(btn_a)
        form.addRow("Molecule A (.xyz)", row_a)

        row_b = QHBoxLayout()
        row_b.addWidget(self.pair_b)
        btn_b = QPushButton("Browse…", tab)
        row_b.addWidget(btn_b)
        form.addRow("Molecule B (.xyz)", row_b)

        form.addRow("Covariance mode", self.pair_cov)
        form.addRow("Regularization", self.pair_reg)
        form.addRow("", self.pair_std)

        run_btn = QPushButton("Run pair comparison", tab)
        form.addRow("", run_btn)

        btn_a.clicked.connect(lambda: self._pick_file(self.pair_a))
        btn_b.clicked.connect(lambda: self._pick_file(self.pair_b))
        run_btn.clicked.connect(self._run_pair)

        self.tabs.addTab(tab, "Pair")

    def _build_library_tab(self):
        tab = QWidget(self)
        form = QFormLayout(tab)

        self.lib_query = QLineEdit(tab)
        self.lib_dir = QLineEdit(tab)
        self.lib_glob = QLineEdit(tab)
        self.lib_glob.setText("*.xyz")
        self.lib_topk = QSpinBox(tab)
        self.lib_topk.setRange(1, 100000)
        self.lib_topk.setValue(10)
        self.lib_cov = QComboBox(tab)
        self.lib_cov.addItems(["full", "diag"])
        self.lib_reg = QDoubleSpinBox(tab)
        self.lib_reg.setRange(1.0e-8, 1.0)
        self.lib_reg.setDecimals(6)
        self.lib_reg.setValue(5.0e-2)
        self.lib_std = QCheckBox("Standardize features", tab)
        self.lib_std.setChecked(True)

        row_q = QHBoxLayout()
        row_q.addWidget(self.lib_query)
        btn_q = QPushButton("Browse…", tab)
        row_q.addWidget(btn_q)
        form.addRow("Query (.xyz)", row_q)

        row_d = QHBoxLayout()
        row_d.addWidget(self.lib_dir)
        btn_d = QPushButton("Browse…", tab)
        row_d.addWidget(btn_d)
        form.addRow("Library directory", row_d)

        form.addRow("Library glob", self.lib_glob)
        form.addRow("Top-k", self.lib_topk)
        form.addRow("Covariance mode", self.lib_cov)
        form.addRow("Regularization", self.lib_reg)
        form.addRow("", self.lib_std)

        run_btn = QPushButton("Run library ranking", tab)
        form.addRow("", run_btn)

        btn_q.clicked.connect(lambda: self._pick_file(self.lib_query))
        btn_d.clicked.connect(lambda: self._pick_directory(self.lib_dir))
        run_btn.clicked.connect(self._run_library)

        self.tabs.addTab(tab, "Library")

    def _pick_file(self, target: QLineEdit):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select XYZ file", str(self.working_dir), "XYZ files (*.xyz);;All files (*)"
        )
        if path:
            target.setText(path)

    def _pick_directory(self, target: QLineEdit):
        path = QFileDialog.getExistingDirectory(
            self, "Select library directory", str(self.working_dir)
        )
        if path:
            target.setText(path)

    def _run_pair(self):
        xyz_a = Path(self.pair_a.text().strip())
        xyz_b = Path(self.pair_b.text().strip())
        if not xyz_a.exists() or not xyz_b.exists():
            QMessageBox.warning(self, "Similarity", "Please select two valid XYZ files.")
            return
        try:
            result = compare_molecules(
                xyz_a,
                xyz_b,
                covariance_mode=self.pair_cov.currentText(),
                regularization=float(self.pair_reg.value()),
                standardize=self.pair_std.isChecked(),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Similarity error", str(exc))
            return

        self._last_result = result
        text = [
            "Pair comparison",
            f"A: {result['xyz_a']}",
            f"B: {result['xyz_b']}",
            f"Similarity exp(-DB): {result['similarity_exp_minus_db']:.8f}",
            f"Bhattacharyya distance: {result['bhattacharyya_distance']:.8f}",
            f"Ring similarity exp(-DB): {result['ring_similarity_exp_minus_db']:.8f}",
            f"Combined similarity: {result['similarity_combined']:.8f}",
            f"Atoms A/B: {result['natoms_a']} / {result['natoms_b']}",
            f"Rings A/B: {result['nrings_a']} / {result['nrings_b']}",
        ]
        top_syn = result.get("explain", {}).get("synthon_feature_terms", [])[:3]
        if top_syn:
            text.append("")
            text.append("Top synthon contributors:")
            for row in top_syn:
                text.append(
                    f" - {row['feature']}: total={row['total_term']:.6f} "
                    f"(mean={row['mean_term']:.6f}, var={row['variance_term']:.6f})"
                )
        self.output.setPlainText("\n".join(text))

    def _run_library(self):
        query = Path(self.lib_query.text().strip())
        lib_dir = Path(self.lib_dir.text().strip())
        if not query.exists() or not lib_dir.exists() or not lib_dir.is_dir():
            QMessageBox.warning(self, "Similarity", "Please set a valid query XYZ and library directory.")
            return
        library = sorted(lib_dir.glob(self.lib_glob.text().strip() or "*.xyz"))
        if not library:
            QMessageBox.warning(self, "Similarity", "No XYZ files found in library with current glob.")
            return
        try:
            result = compare_against_library(
                query,
                library,
                covariance_mode=self.lib_cov.currentText(),
                regularization=float(self.lib_reg.value()),
                standardize=self.lib_std.isChecked(),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Similarity error", str(exc))
            return

        self._last_result = result
        top_k = max(1, int(self.lib_topk.value()))
        ranking = result["ranking"][:top_k]

        lines = [
            "Library ranking",
            f"Query: {result['query_xyz']}",
            f"Compared molecules: {result['library_size']}",
        ]
        skipped = result.get("skipped", [])
        if skipped:
            lines.append(f"Skipped molecules: {len(skipped)}")
        lines.append("")
        for idx, row in enumerate(ranking, start=1):
            lines.append(
                f"{idx:3d}. sim_comb={row['similarity_combined']:.8f}  "
                f"sim_syn={row['similarity_exp_minus_db']:.8f}  "
                f"sim_ring={row['ring_similarity_exp_minus_db']:.8f}  "
                f"{row['xyz']}"
            )
        self.output.setPlainText("\n".join(lines))

    def _save_json_report(self):
        if self._last_result is None:
            QMessageBox.information(self, "Similarity", "No result available yet.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save similarity JSON",
            str(self.working_dir / "similarity.json"),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            Path(path).write_text(json.dumps(self._last_result, indent=2) + "\n", encoding="utf-8")
        except Exception as exc:
            QMessageBox.critical(self, "Save error", str(exc))
            return
        QMessageBox.information(self, "Similarity", f"Saved report:\n{path}")
