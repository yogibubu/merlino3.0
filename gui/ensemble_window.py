from __future__ import annotations

from pathlib import Path
import tomllib

from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from merlino_semiexp import fit_ensemble_job


class EnsembleSEWindow(QDialog):
    def __init__(self, working_dir: Path, parent=None):
        super().__init__(parent)
        self.working_dir = Path(working_dir)
        self._last_result = None
        self._last_outdir: Path | None = None
        self._job_data: dict | None = None

        self.setWindowTitle("MORPHEUS Ensemble SE")
        self.resize(1120, 780)

        root = QVBoxLayout(self)
        title = QLabel("Multi-molecule class-correction refinement from shared coordinate classes")
        title.setStyleSheet("font-weight: 600;")
        root.addWidget(title)

        panel = QWidget(self)
        form = QFormLayout(panel)
        root.addWidget(panel)

        self.job_path = QLineEdit(panel)
        self.out_dir = QLineEdit(panel)
        self.out_dir.setText(str(self.working_dir / "semiexp_ensemble"))

        row_job = QHBoxLayout()
        row_job.addWidget(self.job_path)
        browse_job = QPushButton("Browse…", panel)
        row_job.addWidget(browse_job)
        form.addRow("Ensemble job", row_job)

        row_out = QHBoxLayout()
        row_out.addWidget(self.out_dir)
        browse_out = QPushButton("Browse…", panel)
        row_out.addWidget(browse_out)
        form.addRow("Output directory", row_out)

        actions = QHBoxLayout()
        root.addLayout(actions)
        self.load_btn = QPushButton("Load classes", self)
        self.template_btn = QPushButton("Add CHO templates", self)
        self.auto_prior_btn = QPushButton("Auto priors", self)
        self.add_btn = QPushButton("Add class", self)
        self.remove_btn = QPushButton("Remove class", self)
        self.save_btn = QPushButton("Save edited job", self)
        self.run_btn = QPushButton("Run ensemble fit", self)
        self.compare_btn = QPushButton("Compare priors", self)
        self.open_btn = QPushButton("Open report", self)
        self.open_btn.setEnabled(False)
        self.close_btn = QPushButton("Close", self)
        actions.addWidget(self.load_btn)
        actions.addWidget(self.template_btn)
        actions.addWidget(self.auto_prior_btn)
        actions.addWidget(self.add_btn)
        actions.addWidget(self.remove_btn)
        actions.addWidget(self.save_btn)
        actions.addWidget(self.run_btn)
        actions.addWidget(self.compare_btn)
        actions.addWidget(self.open_btn)
        actions.addStretch()
        actions.addWidget(self.close_btn)

        root.addWidget(QLabel("Class/prior editor", self))
        self.class_table = QTableWidget(0, 11, self)
        self.class_table.setHorizontalHeaderLabels(
            [
                "Name",
                "Kind",
                "Atoms",
                "Patterns",
                "Min value",
                "Max value",
                "Synthons",
                "Zeff",
                "Zeff tol",
                "Prior value",
                "Prior sigma",
            ]
        )
        root.addWidget(self.class_table)

        root.addWidget(QLabel("Fit results", self))
        self.result_table = QTableWidget(0, 7, self)
        self.result_table.setHorizontalHeaderLabels(
            ["Class", "Kind", "Atoms", "Prior", "Correction", "Sigma", "Prior residual"]
        )
        root.addWidget(self.result_table)

        self.output = QTextEdit(self)
        self.output.setReadOnly(True)
        root.addWidget(self.output)

        browse_job.clicked.connect(self._pick_job)
        browse_out.clicked.connect(self._pick_outdir)
        self.load_btn.clicked.connect(self._load_job)
        self.template_btn.clicked.connect(self._add_cho_templates)
        self.auto_prior_btn.clicked.connect(self._apply_auto_priors)
        self.add_btn.clicked.connect(self._add_class)
        self.remove_btn.clicked.connect(self._remove_selected_classes)
        self.save_btn.clicked.connect(self._save_edited_job)
        self.run_btn.clicked.connect(self._run)
        self.compare_btn.clicked.connect(self._compare)
        self.open_btn.clicked.connect(self._open_report)
        self.close_btn.clicked.connect(self.accept)

    def _pick_job(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select ensemble SE job",
            str(self.working_dir),
            "Ensemble jobs (*.toml *.mse-ensemble.toml);;All files (*)",
        )
        if path:
            self.job_path.setText(path)
            if not self.out_dir.text().strip():
                self.out_dir.setText(str(Path(path).parent / "ensemble_out"))
            self._load_job()

    def _pick_outdir(self):
        path = QFileDialog.getExistingDirectory(
            self,
            "Select ensemble output directory",
            self.out_dir.text().strip() or str(self.working_dir),
        )
        if path:
            self.out_dir.setText(path)

    def _run(self):
        job = self._materialize_job_for_run()
        if job is None:
            return
        outdir = Path(self.out_dir.text().strip() or self.working_dir / "semiexp_ensemble")
        try:
            result = fit_ensemble_job(job, outdir=outdir)
        except Exception as exc:
            QMessageBox.critical(self, "MORPHEUS Ensemble SE", str(exc))
            return
        self._last_result = result
        self._last_outdir = outdir
        self._populate_result()
        self.open_btn.setEnabled((outdir / "ensemble_class_corrections.txt").exists())

    def _compare(self):
        from merlino_semiexp import run_ensemble_prior_comparison

        job = self._materialize_job_for_run()
        if job is None:
            return
        outdir = Path(self.out_dir.text().strip() or self.working_dir / "semiexp_ensemble") / "comparison"
        try:
            results = run_ensemble_prior_comparison(job, outdir)
        except Exception as exc:
            QMessageBox.critical(self, "MORPHEUS Ensemble SE", str(exc))
            return
        summary = ["Prior comparison:"]
        for name, result in results.items():
            summary.append(
                f"{name}: classes={len(result.classes)} rank={result.rank} "
                f"acceptance={result.acceptance.status} "
                f"condition={result.condition_number:.6g} "
                f"WRMS={result.weighted_rms_before:.6g}->{result.weighted_rms_after:.6g}"
            )
        self.output.setPlainText("\n".join(summary))
        self._last_outdir = outdir
        self.open_btn.setEnabled(False)

    def _populate_result(self):
        result = self._last_result
        if result is None:
            return
        self.output.setPlainText(result.text)
        self.result_table.setRowCount(len(result.classes))
        for row, item in enumerate(result.classes):
            prior = ""
            if item.prior_value is not None and item.prior_sigma is not None:
                prior = f"{item.prior_value:.6g} +/- {item.prior_sigma:.6g}"
            values = [
                item.name,
                item.kind,
                ",".join(item.atom_symbols),
                prior,
                f"{result.corrections[item.name]:.10g}",
                f"{result.sigma[item.name]:.4g}",
                "" if item.name not in result.prior_residual_after else f"{result.prior_residual_after[item.name]:.6g}",
            ]
            for col, value in enumerate(values):
                self.result_table.setItem(row, col, QTableWidgetItem(value))
        self.result_table.resizeColumnsToContents()

    def _open_report(self):
        if self._last_outdir is None:
            return
        report = self._last_outdir / "ensemble_class_corrections.txt"
        if not report.exists():
            QMessageBox.warning(self, "MORPHEUS Ensemble SE", "Report not found.")
            return
        self.output.setPlainText(report.read_text(encoding="utf-8", errors="replace"))

    def _load_job(self):
        job = Path(self.job_path.text().strip())
        if not job.exists():
            QMessageBox.warning(self, "MORPHEUS Ensemble SE", "Select an existing ensemble job.")
            return
        try:
            self._job_data = tomllib.loads(job.read_text(encoding="utf-8"))
        except Exception as exc:
            QMessageBox.critical(self, "MORPHEUS Ensemble SE", str(exc))
            return
        self._populate_class_editor(self._job_data.get("classes", ()))

    def _populate_class_editor(self, classes):
        self.class_table.setRowCount(0)
        for item in classes or ():
            self._append_class_row(
                str(item.get("name", "")),
                str(item.get("kind", "")),
                _join_items(item.get("atoms", item.get("atom_symbols", ""))),
                _join_items(item.get("patterns", "")),
                "" if item.get("value_min", item.get("min_value")) is None else str(item.get("value_min", item.get("min_value"))),
                "" if item.get("value_max", item.get("max_value")) is None else str(item.get("value_max", item.get("max_value"))),
                _join_items(item.get("synthon_signatures", item.get("synthons", ""))),
                _join_items(item.get("synthon_zeff", item.get("zeff", ""))),
                "" if item.get("synthon_threshold", item.get("zeff_threshold")) is None else str(item.get("synthon_threshold", item.get("zeff_threshold"))),
                "" if item.get("prior_value", item.get("prior")) is None else str(item.get("prior_value", item.get("prior"))),
                "" if item.get("prior_sigma", item.get("sigma")) is None else str(item.get("prior_sigma", item.get("sigma"))),
            )
        self.class_table.resizeColumnsToContents()

    def _append_class_row(
        self,
        name: str = "",
        kind: str = "stretch",
        atoms: str = "",
        patterns: str = "",
        value_min: str = "",
        value_max: str = "",
        synthons: str = "",
        synthon_zeff: str = "",
        synthon_threshold: str = "",
        prior_value: str = "",
        prior_sigma: str = "",
    ):
        row = self.class_table.rowCount()
        self.class_table.insertRow(row)
        for col, value in enumerate(
            (name, kind, atoms, patterns, value_min, value_max, synthons, synthon_zeff, synthon_threshold, prior_value, prior_sigma)
        ):
            self.class_table.setItem(row, col, QTableWidgetItem(value))

    def _add_class(self):
        self._append_class_row()

    def _add_cho_templates(self):
        existing = {_cell(self.class_table, row, 0) for row in range(self.class_table.rowCount())}
        templates = [
            ("CC_short", "stretch", "C,C", "", "", "1.42", "", "", "", "", ""),
            ("CC_long", "stretch", "C,C", "", "1.42", "", "", "", "", "", ""),
            ("CO_carbonyl", "stretch", "C,O", "", "", "1.28", "", "", "", "", ""),
            ("CO_single", "stretch", "C,O", "", "1.28", "", "", "", "", "", ""),
            ("CH_stretch", "stretch", "C,H", "", "", "", "", "", "", "", ""),
            ("CCC_bend", "bend", "C,C,C", "", "", "", "", "", "", "0.0", "1.0e-3"),
            ("CCO_bend", "bend", "C,C,O", "", "", "", "", "", "", "0.0", "1.0e-3"),
            ("COC_bend", "bend", "C,O,C", "", "", "", "", "", "", "0.0", "1.0e-3"),
            ("OCO_bend", "bend", "O,C,O", "", "", "", "", "", "", "0.0", "1.0e-3"),
        ]
        for row in templates:
            if row[0] not in existing:
                self._append_class_row(*row)
        self.class_table.resizeColumnsToContents()

    def _apply_auto_priors(self):
        for row in range(self.class_table.rowCount()):
            kind = _cell(self.class_table, row, 1).lower().replace("-", "_")
            prior_value = _cell(self.class_table, row, 9)
            prior_sigma = _cell(self.class_table, row, 10)
            if prior_value or prior_sigma:
                continue
            sigma = ""
            if kind in {"bend", "angle", "a"}:
                sigma = "1.0e-3"
            elif kind in {"torsion", "dihedral", "d"}:
                sigma = "2.0e-3"
            elif kind in {"out_of_plane", "oop", "u"}:
                sigma = "1.0e-3"
            if sigma:
                self.class_table.setItem(row, 9, QTableWidgetItem("0.0"))
                self.class_table.setItem(row, 10, QTableWidgetItem(sigma))
        self.class_table.resizeColumnsToContents()

    def _remove_selected_classes(self):
        rows = sorted({index.row() for index in self.class_table.selectedIndexes()}, reverse=True)
        for row in rows:
            self.class_table.removeRow(row)

    def _save_edited_job(self):
        job = Path(self.job_path.text().strip())
        if not job:
            return
        if self._job_data is None:
            self._load_job()
        if self._job_data is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save edited ensemble job",
            str(job),
            "Ensemble jobs (*.toml *.mse-ensemble.toml);;All files (*)",
        )
        if not path:
            return
        target = Path(path)
        target.write_text(
            _ensemble_job_text(self._job_data, self._classes_from_editor(), source_path=job),
            encoding="utf-8",
        )
        self.job_path.setText(str(target))
        self._job_data = tomllib.loads(target.read_text(encoding="utf-8"))

    def _materialize_job_for_run(self) -> Path | None:
        job = Path(self.job_path.text().strip())
        if not job.exists():
            QMessageBox.warning(self, "MORPHEUS Ensemble SE", "Select an existing ensemble job.")
            return None
        if self._job_data is None:
            self._load_job()
        if self._job_data is None:
            return None
        outdir = Path(self.out_dir.text().strip() or self.working_dir / "semiexp_ensemble")
        outdir.mkdir(parents=True, exist_ok=True)
        target = outdir / "edited_ensemble_job.mse-ensemble.toml"
        target.write_text(
            _ensemble_job_text(self._job_data, self._classes_from_editor(), source_path=job),
            encoding="utf-8",
        )
        return target

    def _classes_from_editor(self) -> list[dict]:
        classes: list[dict] = []
        for row in range(self.class_table.rowCount()):
            name = _cell(self.class_table, row, 0)
            kind = _cell(self.class_table, row, 1)
            atoms = _split_items(_cell(self.class_table, row, 2))
            patterns = _split_items(_cell(self.class_table, row, 3))
            value_min = _cell(self.class_table, row, 4)
            value_max = _cell(self.class_table, row, 5)
            synthons = _split_items(_cell(self.class_table, row, 6))
            synthon_zeff = _split_items(_cell(self.class_table, row, 7))
            synthon_threshold = _cell(self.class_table, row, 8)
            prior_value = _cell(self.class_table, row, 9)
            prior_sigma = _cell(self.class_table, row, 10)
            if not name:
                continue
            item = {"name": name}
            if kind:
                item["kind"] = kind
            if atoms:
                item["atoms"] = atoms
            if patterns:
                item["patterns"] = patterns
            if value_min:
                item["value_min"] = float(value_min)
            if value_max:
                item["value_max"] = float(value_max)
            if synthons:
                item["synthon_signatures"] = synthons
            if synthon_zeff:
                item["synthon_zeff"] = [float(value) for value in synthon_zeff]
            if synthon_threshold:
                item["synthon_threshold"] = float(synthon_threshold)
            if prior_value:
                item["prior_value"] = float(prior_value)
            if prior_sigma:
                item["prior_sigma"] = float(prior_sigma)
            classes.append(item)
        return classes


def _cell(table: QTableWidget, row: int, col: int) -> str:
    item = table.item(row, col)
    return "" if item is None else item.text().strip()


def _split_items(text: str) -> list[str]:
    return [part.strip() for part in text.replace("|", ",").replace(";", ",").split(",") if part.strip()]


def _join_items(value) -> str:
    if isinstance(value, str):
        return value
    try:
        return ",".join(str(item) for item in value)
    except TypeError:
        return ""


def _quote(value: str) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def _toml_array(values: list[str]) -> str:
    return "[" + ", ".join(_quote(value) for value in values) + "]"


def _ensemble_job_text(data: dict, classes: list[dict], *, source_path: Path | None = None) -> str:
    lines = [
        'schema = "merlino.semiexp.ensemble.v1"',
        f"title = {_quote(str(data.get('title', 'ensemble class-correction fit')))}",
        "",
        "[fit]",
    ]
    fit = data.get("fit", {}) if isinstance(data.get("fit", {}), dict) else {}
    lines.append(f"step = {float(fit.get('step', 1.0e-4)):.12g}")
    lines.append(f"rcond = {float(fit.get('rcond', 1.0e-10)):.12g}")
    acceptance = data.get("acceptance", {}) if isinstance(data.get("acceptance", {}), dict) else {}
    lines.extend(["", "[acceptance]"])
    lines.append(f"require_full_rank = {_toml_bool(acceptance.get('require_full_rank', True))}")
    lines.append(f"max_condition_number = {float(acceptance.get('max_condition_number', 1.0e8)):.12g}")
    lines.append(
        "min_residual_degrees_of_freedom = "
        f"{int(acceptance.get('min_residual_degrees_of_freedom', 1))}"
    )
    lines.append(f"min_molecule_support = {int(acceptance.get('min_molecule_support', 2))}")
    lines.append(
        "high_correlation_review_threshold = "
        f"{float(acceptance.get('high_correlation_review_threshold', 0.98)):.12g}"
    )
    lines.append(
        "high_correlation_reject_threshold = "
        f"{float(acceptance.get('high_correlation_reject_threshold', 0.9999)):.12g}"
    )
    for molecule in data.get("molecules", ()):
        lines.extend(["", "[[molecules]]"])
        lines.append(f"name = {_quote(str(molecule.get('name', molecule.get('label', 'molecule'))))}")
        molecule_job = str(molecule.get("job", ""))
        if source_path is not None and molecule_job and not Path(molecule_job).is_absolute():
            molecule_job = str((source_path.parent / molecule_job).resolve())
        lines.append(f"job = {_quote(molecule_job)}")
    for item in classes:
        lines.extend(["", "[[classes]]"])
        lines.append(f"name = {_quote(str(item['name']))}")
        if item.get("kind"):
            lines.append(f"kind = {_quote(str(item['kind']))}")
        if item.get("atoms"):
            lines.append(f"atoms = {_toml_array(list(item['atoms']))}")
        if item.get("patterns"):
            lines.append(f"patterns = {_toml_array(list(item['patterns']))}")
        if "value_min" in item:
            lines.append(f"value_min = {float(item['value_min']):.12g}")
        if "value_max" in item:
            lines.append(f"value_max = {float(item['value_max']):.12g}")
        if item.get("synthon_signatures"):
            lines.append(f"synthon_signatures = {_toml_array(list(item['synthon_signatures']))}")
        if item.get("synthon_zeff"):
            lines.append("synthon_zeff = [" + ", ".join(f"{float(value):.12g}" for value in item["synthon_zeff"]) + "]")
        if "synthon_threshold" in item:
            lines.append(f"synthon_threshold = {float(item['synthon_threshold']):.12g}")
        if "prior_value" in item:
            lines.append(f"prior_value = {float(item['prior_value']):.12g}")
        if "prior_sigma" in item:
            lines.append(f"prior_sigma = {float(item['prior_sigma']):.12g}")
    return "\n".join(lines) + "\n"


def _toml_bool(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    text = str(value).strip().lower()
    return "true" if text in {"1", "true", "yes", "y", "on"} else "false"
