from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QCheckBox,
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
from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices

from merlino_fit.survibfit.fragment_pipeline import run_fragment_pipeline
from merlino_fit.survibfit.fragment_delta_correction import (
    prepare_delta_workflow,
    prepare_hpcs2_delta_workflow,
    prepare_low_level_gaussian_inputs,
    prepare_high_level_curation_jobs,
    apply_delta_correction,
)


class FragmentPipelineWindow(QDialog):
    def __init__(self, working_dir: Path, parent=None):
        super().__init__(parent)
        self.working_dir = Path(working_dir)
        self._last_report = None

        self.setWindowTitle("Fragment Pipeline (MVP)")
        self.resize(860, 620)

        root = QVBoxLayout(self)
        title = QLabel(
            "Fragment similarity + GAP detection (SE preferred over PCS2 on duplicate molecules; "
            "PCS2 target uses HPCS2 base geometry)"
        )
        title.setStyleSheet("font-weight: 600;")
        root.addWidget(title)

        panel = QWidget(self)
        form = QFormLayout(panel)
        root.addWidget(panel)

        self.query_xyz = QLineEdit(panel)
        self.se_dir = QLineEdit(panel)
        self.pcs2_dir = QLineEdit(panel)
        self.hpcs2_dir = QLineEdit(panel)
        self.out_dir = QLineEdit(panel)
        self.top_k = QSpinBox(panel)
        self.top_k.setRange(1, 200)
        self.top_k.setValue(5)
        self.max_ring_size_delta = QSpinBox(panel)
        self.max_ring_size_delta.setRange(0, 12)
        self.max_ring_size_delta.setValue(2)
        self.chemical_min_score = QDoubleSpinBox(panel)
        self.chemical_min_score.setRange(0.0, 1.0)
        self.chemical_min_score.setDecimals(3)
        self.chemical_min_score.setValue(0.15)
        self.reuse_penalty = QDoubleSpinBox(panel)
        self.reuse_penalty.setRange(0.0, 1.0)
        self.reuse_penalty.setDecimals(3)
        self.reuse_penalty.setValue(0.04)
        self.low_score_threshold = QDoubleSpinBox(panel)
        self.low_score_threshold.setRange(0.0, 1.0)
        self.low_score_threshold.setDecimals(3)
        self.low_score_threshold.setValue(0.45)
        self.gap_threshold = QDoubleSpinBox(panel)
        self.gap_threshold.setRange(0.0, 1.0)
        self.gap_threshold.setDecimals(3)
        self.gap_threshold.setValue(0.75)
        self.low_level_route = QLineEdit(panel)
        self.low_level_route.setText("#p B3LYP/6-31G(d) Opt")
        self.low_level_charge = QSpinBox(panel)
        self.low_level_charge.setRange(-10, 10)
        self.low_level_charge.setValue(0)
        self.low_level_mult = QSpinBox(panel)
        self.low_level_mult.setRange(1, 10)
        self.low_level_mult.setValue(1)
        self.low_level_nproc = QSpinBox(panel)
        self.low_level_nproc.setRange(1, 64)
        self.low_level_nproc.setValue(8)
        self.low_level_mem = QLineEdit(panel)
        self.low_level_mem.setText("8GB")
        self.high_level_route = QLineEdit(panel)
        self.high_level_route.setText("#p wB97XD/def2TZVP Opt")
        self.high_level_nproc = QSpinBox(panel)
        self.high_level_nproc.setRange(1, 64)
        self.high_level_nproc.setValue(8)
        self.high_level_mem = QLineEdit(panel)
        self.high_level_mem.setText("16GB")
        self.target_pcs2 = QCheckBox("Target level: PCS2", panel)
        self.target_pcs2.setToolTip(
            "PCS2 target: match against PCS2 fragment library and correct HPCS2 base geometry."
        )
        self.target_note = QLabel(
            "Note: PCS2 uses fragment library (high-level); HPCS2 is full-molecule base (low-level).",
            panel,
        )
        self.target_note.setStyleSheet("color: #555;")

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

        row_hpcs2 = QHBoxLayout()
        row_hpcs2.addWidget(self.hpcs2_dir)
        btn_hpcs2 = QPushButton("Browse…", panel)
        row_hpcs2.addWidget(btn_hpcs2)
        form.addRow("HPCS2 library dir", row_hpcs2)

        row_out = QHBoxLayout()
        row_out.addWidget(self.out_dir)
        btn_out = QPushButton("Browse…", panel)
        row_out.addWidget(btn_out)
        form.addRow("Output dir", row_out)

        form.addRow("Top-k per fragment", self.top_k)
        form.addRow("Max ring-size delta", self.max_ring_size_delta)
        form.addRow("Chemical min score", self.chemical_min_score)
        form.addRow("Reuse penalty", self.reuse_penalty)
        form.addRow("Low-score threshold", self.low_score_threshold)
        form.addRow("GAP threshold", self.gap_threshold)
        form.addRow("LL Gaussian route", self.low_level_route)
        form.addRow("LL charge", self.low_level_charge)
        form.addRow("LL multiplicity", self.low_level_mult)
        form.addRow("LL nproc", self.low_level_nproc)
        form.addRow("LL memory", self.low_level_mem)
        form.addRow("HL Gaussian route", self.high_level_route)
        form.addRow("HL nproc", self.high_level_nproc)
        form.addRow("HL memory", self.high_level_mem)
        form.addRow("", self.target_pcs2)
        form.addRow("", self.target_note)

        run_btn = QPushButton("Run fragment pipeline", panel)
        prep_ll_btn = QPushButton("Prepare LL fragment Gaussian inputs", panel)
        prep_hl_btn = QPushButton("Export HL curation queue jobs", panel)
        apply_hpcs2_btn = QPushButton("Apply PCS2 delta correction", panel)
        form.addRow("", run_btn)
        form.addRow("", prep_ll_btn)
        form.addRow("", prep_hl_btn)
        form.addRow("", apply_hpcs2_btn)

        self.output = QTextEdit(self)
        self.output.setReadOnly(True)
        root.addWidget(self.output)

        footer = QHBoxLayout()
        root.addLayout(footer)
        save_btn = QPushButton("Save JSON report")
        save_queue_btn = QPushButton("Save HL curation queue")
        close_btn = QPushButton("Close")
        footer.addWidget(save_btn)
        footer.addWidget(save_queue_btn)
        footer.addStretch()
        footer.addWidget(close_btn)

        btn_q.clicked.connect(lambda: self._pick_file(self.query_xyz))
        btn_se.clicked.connect(lambda: self._pick_dir(self.se_dir))
        btn_pcs2.clicked.connect(lambda: self._pick_dir(self.pcs2_dir))
        btn_hpcs2.clicked.connect(lambda: self._pick_dir(self.hpcs2_dir))
        btn_out.clicked.connect(lambda: self._pick_dir(self.out_dir))
        run_btn.clicked.connect(self._run_pipeline)
        prep_ll_btn.clicked.connect(self._prepare_low_level_inputs)
        prep_hl_btn.clicked.connect(self._prepare_high_level_queue_jobs)
        apply_hpcs2_btn.clicked.connect(self._apply_pcs2_correction)
        save_btn.clicked.connect(self._save_json_report)
        save_queue_btn.clicked.connect(self._save_curation_queue)
        close_btn.clicked.connect(self.accept)

        # Defaults: prefer project-local libraries, fallback to parent-level layout.
        se_local = self.working_dir / "projects" / "se_library"
        pcs2_local = self.working_dir / "projects" / "pcs2_library"
        hpcs2_local = self.working_dir / "projects" / "hpcs2_library"
        se_parent = self.working_dir.parent / "projects" / "se_library"
        pcs2_parent = self.working_dir.parent / "projects" / "pcs2_library"
        hpcs2_parent = self.working_dir.parent / "projects" / "hpcs2_library"
        self.se_dir.setText(str(se_local if se_local.is_dir() else se_parent))
        self.pcs2_dir.setText(str(pcs2_local if pcs2_local.is_dir() else pcs2_parent))
        self.hpcs2_dir.setText(str(hpcs2_local if hpcs2_local.is_dir() else hpcs2_parent))
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
        if self.target_pcs2.isChecked():
            if not pcs2_dir.is_dir():
                QMessageBox.warning(
                    self, "Fragment pipeline", "Please select a valid PCS2 library directory."
                )
                return
        elif not se_dir.is_dir() or not pcs2_dir.is_dir():
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
                max_ring_size_delta=int(self.max_ring_size_delta.value()),
                chemical_min_score=float(self.chemical_min_score.value()),
                reuse_penalty=float(self.reuse_penalty.value()),
                low_score_threshold=float(self.low_score_threshold.value()),
                use_pcs2_only=self.target_pcs2.isChecked(),
                library_label="PCS2/HPCS2" if self.target_pcs2.isChecked() else "SE/PCS2",
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
            f"Max ring-size delta: {report.get('max_ring_size_delta', 2)}",
            f"Chemical min score: {report.get('chemical_min_score', 0.25):.3f}",
            f"Reuse penalty: {report.get('reuse_penalty', 0.08):.3f}",
            f"Low-score threshold: {report.get('low_score_threshold', 0.45):.3f}",
            f"Low-score fragments: {report.get('low_score_fragments', 0)}",
            f"GAP fragments: {len(report['to_curate'])}",
            "",
            "Top fragment statuses:",
        ]
        for row in report["fragments"][:10]:
            msg = (
                f" - {row['query_fragment_id']}: best={row['best_similarity']:.6f} "
                f"status={row['status']}"
            )
            if row.get("ranking"):
                r0 = row["ranking"][0]
                best_label = r0.get("candidate_fragment_label", r0.get("candidate_fragment_id", ""))
                msg += (
                    f" -> {best_label}"
                    f" (feat={r0.get('similarity_feature', 0.0):.6f}, "
                    f"3D={r0.get('similarity_3d', 0.0):.6f}, "
                    f"chem={r0.get('chemical_compatibility', 0.0):.6f}, "
                    f"size_d={r0.get('ring_size_delta', 0)}, "
                    f"rmsd={r0.get('ring_rmsd_3d', 0.0):.4f})"
                )
            if row.get("is_low_score"):
                sug = row.get("suggested_new_fragment") or {}
                kind = sug.get("query_fragment_kind", "fragment")
                size = sug.get("query_core_size", "?")
                ref = (sug.get("closest_library_reference") or {}).get("candidate_fragment_label", "-")
                msg += (
                    f" [suggest NEW HL fragment: query-derived {kind}_{size}, ref={ref}, "
                    f"rel={row.get('best_reliability', 0.0):.2f}]"
                )
            lines.append(msg)
        lines += [
            "",
            f"Saved: {out_dir / 'fragment_pipeline.json'}",
            f"Saved: {out_dir / 'to_curate.json'}",
            f"Saved: {out_dir / 'fragment_pipeline.md'}",
        ]
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
        assembly = report.get("global_assembly", {})
        lines += ["", "Global assembly:"]
        if assembly.get("available"):
            lines.append(
                f" score={assembly.get('score', 0.0):.6f} "
                f"ind_mean={assembly.get('individual_mean', 0.0):.6f} "
                f"pair={assembly.get('pairwise_mean', 0.0):.6f} "
                f"reuse={assembly.get('reuse_ratio', 0.0):.6f}"
            )
            for a in assembly.get("assignment", [])[:10]:
                lines.append(
                    f" - {a['query_fragment_id']} -> {a['candidate_fragment_label']} "
                    f"(sim={a['similarity']:.6f}, chem={a['chemical_compatibility']:.6f})"
                )
        else:
            lines.append(f" not available: {assembly.get('reason', 'unknown')}")
        queue = report.get("high_level_curation_queue", [])
        lines += ["", "HL curation queue:"]
        if queue:
            for q in queue[:20]:
                ref = q.get("closest_library_reference", {})
                lines.append(
                    f" - {q.get('query_fragment_id')} -> new {q.get('query_fragment_kind')}_{q.get('query_core_size')} "
                    f"(best={q.get('best_similarity', 0.0):.6f}, ref={ref.get('candidate_fragment_label', '-')})"
                )
        else:
            lines.append(" - none")
        self.output.setPlainText("\n".join(lines))

    def _apply_pcs2_correction(self):
        query = Path(self.query_xyz.text().strip())
        out_dir = Path(self.out_dir.text().strip())
        hpcs2_dir = Path(self.hpcs2_dir.text().strip())
        report_json = out_dir / "fragment_pipeline.json"
        if not query.exists():
            QMessageBox.warning(self, "PCS2 correction", "Please select a valid query XYZ file.")
            return
        if not hpcs2_dir.is_dir():
            QMessageBox.warning(self, "PCS2 correction", "Please select a valid HPCS2 library directory.")
            return
        if not report_json.exists():
            QMessageBox.warning(
                self,
                "PCS2 correction",
                "Run fragment pipeline first so fragment_pipeline.json is available.",
            )
            return

        delta_dir = out_dir / "delta_bundle_hpcs2"
        try:
            manifest = prepare_hpcs2_delta_workflow(query, report_json, hpcs2_dir, delta_dir)
        except Exception as exc:
            QMessageBox.critical(self, "PCS2 correction error", str(exc))
            return

        if not manifest.get("entries"):
            QMessageBox.warning(
                self,
                "PCS2 correction",
                "No valid PCS2->HPCS2 fragment pairs found. Check libraries and rerun.",
            )
            return

        out_xyz = self.working_dir / f"{query.stem}.pcs2.xyz"
        try:
            meta = apply_delta_correction(query, delta_dir / "delta_manifest.json", out_xyz)
        except Exception as exc:
            QMessageBox.critical(self, "PCS2 correction error", str(exc))
            return
        try:
            from merlino_fit.survibfit.fragment_pipeline import write_fragment_view_html

            report = json.loads(report_json.read_text(encoding="utf-8"))
            write_fragment_view_html(
                query,
                report.get("fragments", []),
                out_dir,
                corrected_xyz=out_xyz,
                corrected_label="PCS2",
            )
        except Exception:
            pass

        lines = self.output.toPlainText().splitlines() if self.output.toPlainText().strip() else []
        lines += [
            "",
            "PCS2 delta correction completed",
            f"Delta manifest: {delta_dir / 'delta_manifest.json'}",
            f"Fragments used: {meta.get('entries_used', 0)} / {meta.get('entries_total', 0)}",
            f"Corrected XYZ: {out_xyz}",
        ]
        skipped = manifest.get("skipped", [])
        if skipped:
            lines.append(f"Skipped fragments: {len(skipped)}")
        self.output.setPlainText("\n".join(lines))

    def _prepare_low_level_inputs(self):
        query = Path(self.query_xyz.text().strip())
        out_dir = Path(self.out_dir.text().strip())
        report_json = out_dir / "fragment_pipeline.json"
        if not query.exists():
            QMessageBox.warning(self, "Fragment pipeline", "Please select a valid query XYZ file.")
            return
        if not report_json.exists():
            QMessageBox.warning(
                self,
                "Fragment pipeline",
                "Run fragment pipeline first so fragment_pipeline.json is available.",
            )
            return
        route = self.low_level_route.text().strip()
        if not route:
            QMessageBox.warning(self, "Fragment pipeline", "Please provide a low-level Gaussian route.")
            return

        delta_dir = out_dir / "delta_bundle"
        try:
            prep = prepare_delta_workflow(query, report_json, delta_dir)
            ll = prepare_low_level_gaussian_inputs(
                delta_dir / "delta_manifest.json",
                route=route,
                charge=int(self.low_level_charge.value()),
                multiplicity=int(self.low_level_mult.value()),
                nproc=int(self.low_level_nproc.value()),
                mem=self.low_level_mem.text().strip() or "8GB",
            )
        except Exception as exc:
            QMessageBox.critical(self, "Low-level prep error", str(exc))
            return

        lines = self.output.toPlainText().splitlines() if self.output.toPlainText().strip() else []
        lines += [
            "",
            "Low-level Gaussian prep completed",
            f"Delta manifest: {delta_dir / 'delta_manifest.json'}",
            f"Fragments prepared: {len(prep.get('entries', []))}",
            f"Gaussian inputs generated: {ll.get('generated_inputs', 0)}",
            f"Route: {route}",
        ]
        self.output.setPlainText("\n".join(lines))

    def _prepare_high_level_queue_jobs(self):
        query = Path(self.query_xyz.text().strip())
        out_dir = Path(self.out_dir.text().strip())
        report_json = out_dir / "fragment_pipeline.json"
        if not query.exists():
            QMessageBox.warning(self, "Fragment pipeline", "Please select a valid query XYZ file.")
            return
        if not report_json.exists():
            QMessageBox.warning(
                self,
                "Fragment pipeline",
                "Run fragment pipeline first so fragment_pipeline.json is available.",
            )
            return
        route = self.high_level_route.text().strip()
        if not route:
            QMessageBox.warning(self, "Fragment pipeline", "Please provide a high-level Gaussian route.")
            return

        delta_dir = out_dir / "delta_bundle"
        hl_jobs_dir = out_dir / "high_level_curation_jobs"
        try:
            prepare_delta_workflow(query, report_json, delta_dir)
            out = prepare_high_level_curation_jobs(
                delta_dir / "delta_manifest.json",
                hl_jobs_dir,
                route=route,
                charge=int(self.low_level_charge.value()),
                multiplicity=int(self.low_level_mult.value()),
                nproc=int(self.high_level_nproc.value()),
                mem=self.high_level_mem.text().strip() or "16GB",
            )
        except Exception as exc:
            QMessageBox.critical(self, "HL queue prep error", str(exc))
            return

        lines = self.output.toPlainText().splitlines() if self.output.toPlainText().strip() else []
        lines += [
            "",
            "High-level curation queue export completed",
            f"Jobs generated: {len(out.get('jobs', []))}",
            f"Jobs dir: {hl_jobs_dir}",
            f"Queue JSON: {out.get('queue_path')}",
            f"Route: {route}",
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

    def _save_curation_queue(self):
        if self._last_report is None:
            QMessageBox.information(self, "Fragment pipeline", "No result available yet.")
            return
        queue = self._last_report.get("high_level_curation_queue", [])
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save high-level curation queue",
            str(self.working_dir / "high_level_curation_queue.json"),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            Path(path).write_text(json.dumps(queue, indent=2) + "\n", encoding="utf-8")
        except Exception as exc:
            QMessageBox.critical(self, "Save error", str(exc))
            return
        QMessageBox.information(self, "Fragment pipeline", f"Saved queue:\n{path}")
