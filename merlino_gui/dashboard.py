from __future__ import annotations

import json
from pathlib import Path
import shlex

from PySide6.QtCore import QProcess, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QTextCursor
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QPushButton,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from merlino_semiexp import (
    ParameterClassConstraint,
    QMParameterPredicate,
    SemiexperimentalFitRequest,
    preview_semiexperimental_conditioning,
    preview_semiexperimental_gics,
    read_observations,
    validate_semiexperimental_request,
)

from .manifest_browser import ManifestBrowserWindow
from .workflow_registry import WorkflowSpec, default_workflows


class DashboardWindow(QMainWindow):
    """Thin Merlino4 dashboard for workflow-oriented navigation."""

    def __init__(
        self,
        workdir: Path,
        workflows: list[WorkflowSpec] | None = None,
        parent=None,
    ):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.workflows = workflows or default_workflows()
        self.selected_backends = {
            workflow.workflow_id: workflow.default_backend for workflow in self.workflows
        }
        self._semiexp_process: QProcess | None = None

        self.setWindowTitle("Merlino 4.0")
        self.resize(1080, 720)
        self._build_menus()
        self._build_ui()

    def _build_menus(self) -> None:
        menubar = self.menuBar()
        menubar.clear()
        project_menu = menubar.addMenu("Project")
        project_menu.addAction("Manifest Browser", self.open_manifest_browser)
        project_menu.addSeparator()
        project_menu.addAction("Exit", self.close)

        grouped = _group_workflows(self.workflows)
        for category in ("Structure", "Coordinates", "Vibrations", "Dynamics", "Project"):
            if category not in grouped:
                continue
            menu: QMenu = menubar.addMenu(category)
            for workflow in grouped[category]:
                action = menu.addAction(workflow.title)
                action.triggered.connect(lambda checked=False, wid=workflow.workflow_id: self.select_workflow(wid))

    def _build_ui(self) -> None:
        central = QWidget(self)
        layout = QVBoxLayout(central)
        self.setCentralWidget(central)

        header = QLabel("Merlino 4.0 Workflow Dashboard")
        header.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(header)

        self.status_label = QLabel(f"Project workdir: {self.workdir}")
        self.status_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.status_label)

        manifest_button = QPushButton("Open Manifest Browser")
        manifest_button.clicked.connect(self.open_manifest_browser)
        layout.addWidget(manifest_button)

        backend_row = QHBoxLayout()
        backend_row.addWidget(QLabel("Backend for selected step:"))
        self.backend_selector = QComboBox()
        self.backend_selector.currentTextChanged.connect(self._backend_changed)
        backend_row.addWidget(self.backend_selector, stretch=1)
        layout.addLayout(backend_row)

        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter, stretch=1)

        self.workflow_list = QTreeWidget()
        self.workflow_list.setHeaderLabels(["Workflow", "Service"])
        splitter.addWidget(self.workflow_list)

        self.detail_view = QTextEdit()
        self.detail_view.setReadOnly(True)
        splitter.addWidget(self.detail_view)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)

        self.semiexp_panel = self._build_semiexp_panel()
        layout.addWidget(self.semiexp_panel)
        self.semiexp_panel.hide()

        for category, workflows in _group_workflows(self.workflows).items():
            parent = QTreeWidgetItem([category, ""])
            parent.setFlags(parent.flags() & ~Qt.ItemIsSelectable)
            self.workflow_list.addTopLevelItem(parent)
            for workflow in workflows:
                item = QTreeWidgetItem([workflow.title, workflow.service])
                item.setData(0, Qt.UserRole, workflow.workflow_id)
                parent.addChild(item)
            parent.setExpanded(True)

        self.workflow_list.currentItemChanged.connect(self._show_workflow)
        self.workflow_list.resizeColumnToContents(0)
        if self.workflows:
            self.select_workflow(self.workflows[0].workflow_id)

    def open_manifest_browser(self) -> None:
        self.manifest_browser = ManifestBrowserWindow(self.workdir, parent=self)
        self.manifest_browser.show()
        self.manifest_browser.raise_()
        self.manifest_browser.activateWindow()

    def select_workflow(self, workflow_id: str) -> None:
        for top_idx in range(self.workflow_list.topLevelItemCount()):
            parent = self.workflow_list.topLevelItem(top_idx)
            for child_idx in range(parent.childCount()):
                child = parent.child(child_idx)
                if child.data(0, Qt.UserRole) == workflow_id:
                    self.workflow_list.setCurrentItem(child)
                    return

    def _show_workflow(self, current: QTreeWidgetItem | None, previous=None) -> None:
        if current is None:
            self.detail_view.clear()
            return
        workflow_id = current.data(0, Qt.UserRole)
        if workflow_id is None:
            return
        workflow = next(item for item in self.workflows if item.workflow_id == workflow_id)
        self._sync_backend_selector(workflow)
        backend = self.selected_backends.get(workflow.workflow_id, workflow.default_backend)
        self.detail_view.setPlainText(workflow_detail_text(workflow, selected_backend=backend, workdir=self.workdir))
        self.semiexp_panel.setVisible(workflow.workflow_id == "semiexp_geometry")
        if workflow.workflow_id == "semiexp_geometry":
            self._update_semiexp_preview()

    def _sync_backend_selector(self, workflow: WorkflowSpec) -> None:
        self.backend_selector.blockSignals(True)
        self.backend_selector.clear()
        self.backend_selector.addItems(workflow.backends)
        selected = self.selected_backends.get(workflow.workflow_id, workflow.default_backend)
        index = self.backend_selector.findText(selected)
        self.backend_selector.setCurrentIndex(index if index >= 0 else 0)
        self.backend_selector.setEnabled(len(workflow.backends) > 1)
        self.backend_selector.blockSignals(False)

    def _backend_changed(self, backend: str) -> None:
        current = self.workflow_list.currentItem()
        if current is None:
            return
        workflow_id = current.data(0, Qt.UserRole)
        if workflow_id is None:
            return
        self.selected_backends[str(workflow_id)] = backend
        workflow = next(item for item in self.workflows if item.workflow_id == workflow_id)
        self.detail_view.setPlainText(workflow_detail_text(workflow, selected_backend=backend, workdir=self.workdir))
        if workflow.workflow_id == "semiexp_geometry":
            self._update_semiexp_preview()

    def _build_semiexp_panel(self) -> QGroupBox:
        panel = QGroupBox("Semiexperimental Geometry Run")
        layout = QVBoxLayout(panel)
        tabs = QTabWidget()
        layout.addWidget(tabs)

        input_tab = QWidget()
        input_layout = QVBoxLayout(input_tab)
        form = QFormLayout()
        input_layout.addLayout(form)

        self.semiexp_xyz = _path_row(self, form, "Parent XYZ:", "Select parent XYZ", "*.xyz")
        self.semiexp_observations = _path_row(self, form, "Isotopologues:", "Select isotopologue observations", "*.toml *.json *.csv")
        self.semiexp_outdir = _directory_row(self, form, "Output directory:")

        self.semiexp_iso_table = QTableWidget(0, 14)
        self.semiexp_iso_table.setHorizontalHeaderLabels([
            "label",
            "substitutions",
            "A_MHz",
            "B_MHz",
            "C_MHz",
            "dvib_A",
            "dvib_B",
            "dvib_C",
            "delec_A",
            "delec_B",
            "delec_C",
            "sigma_A",
            "sigma_B",
            "sigma_C",
        ])
        self.semiexp_iso_table.setMaximumHeight(150)
        input_layout.addWidget(QLabel("Isotopologue editor (optional, writes TOML):"))
        input_layout.addWidget(self.semiexp_iso_table)

        option_tab = QWidget()
        option_form = QFormLayout(option_tab)
        self.semiexp_observable = QComboBox()
        self.semiexp_observable.addItems(["moments", "rotational_constants", "auto"])
        self.semiexp_observable.currentTextChanged.connect(lambda _text: self._update_semiexp_preview())
        option_form.addRow("Fit target:", self.semiexp_observable)

        self.semiexp_components = QComboBox()
        self.semiexp_components.addItems(["auto", "ABC", "AB", "AC", "BC"])
        self.semiexp_components.currentTextChanged.connect(lambda _text: self._update_semiexp_preview())
        option_form.addRow("Rotational components:", self.semiexp_components)

        self.semiexp_fixed = QLineEdit()
        self.semiexp_fixed.setPlaceholderText("e.g. GIC001, angle(2,1,3)")
        self.semiexp_fixed.textChanged.connect(lambda _text: self._update_semiexp_preview())
        option_form.addRow("Fixed GIC patterns:", self.semiexp_fixed)

        self.semiexp_qm = QLineEdit()
        self.semiexp_qm.setPlaceholderText("pattern:value:sigma[:source]; repeat with semicolons")
        self.semiexp_qm.textChanged.connect(lambda _text: self._update_semiexp_preview())
        option_form.addRow("QM predicates:", self.semiexp_qm)

        self.semiexp_classes = QLineEdit()
        self.semiexp_classes.setPlaceholderText("CH:shared:bond(1,2)|bond(1,3); XYH:fixed:angle")
        self.semiexp_classes.textChanged.connect(lambda _text: self._update_semiexp_preview())
        option_form.addRow("Parameter classes:", self.semiexp_classes)

        self.semiexp_prune_condition = QDoubleSpinBox()
        self.semiexp_prune_condition.setRange(0.0, 1.0e9)
        self.semiexp_prune_condition.setDecimals(1)
        self.semiexp_prune_condition.setValue(200.0)
        self.semiexp_prune_condition.setSpecialValueText("disabled")
        self.semiexp_prune_condition.valueChanged.connect(lambda _value: self._update_semiexp_preview())
        option_form.addRow("Prune condition target:", self.semiexp_prune_condition)

        preview_tab = QWidget()
        preview_layout = QVBoxLayout(preview_tab)
        self.semiexp_preview_table = QTableWidget(0, 5)
        self.semiexp_preview_table.setHorizontalHeaderLabels(["label", "type", "atoms", "suggested_class", "state"])
        preview_layout.addWidget(self.semiexp_preview_table)

        self.semiexp_command = QTextEdit()
        self.semiexp_command.setReadOnly(True)
        self.semiexp_command.setMaximumHeight(90)
        run_tab = QWidget()
        run_layout = QVBoxLayout(run_tab)
        run_layout.addWidget(QLabel("Command preview:"))
        run_layout.addWidget(self.semiexp_command)

        buttons = QHBoxLayout()
        add_iso_button = QPushButton("Add Isotopologue")
        add_iso_button.clicked.connect(self.add_semiexp_isotopologue_row)
        buttons.addWidget(add_iso_button)
        save_iso_button = QPushButton("Save TOML")
        save_iso_button.clicked.connect(self.save_semiexp_observations_toml)
        buttons.addWidget(save_iso_button)
        preview_gic_button = QPushButton("Preview GIC")
        preview_gic_button.clicked.connect(self.preview_semiexp_gics)
        buttons.addWidget(preview_gic_button)
        condition_button = QPushButton("Check Conditioning")
        condition_button.clicked.connect(self.preview_semiexp_conditioning)
        buttons.addWidget(condition_button)
        validate_button = QPushButton("Validate Input")
        validate_button.clicked.connect(self.validate_semiexp_input)
        buttons.addWidget(validate_button)
        suggest_classes_button = QPushButton("Suggest Classes")
        suggest_classes_button.clicked.connect(self.suggest_semiexp_classes)
        buttons.addWidget(suggest_classes_button)
        save_preset_button = QPushButton("Save Preset")
        save_preset_button.clicked.connect(self.save_semiexp_preset)
        buttons.addWidget(save_preset_button)
        load_preset_button = QPushButton("Load Preset")
        load_preset_button.clicked.connect(self.load_semiexp_preset)
        buttons.addWidget(load_preset_button)
        open_report_button = QPushButton("Open Report")
        open_report_button.clicked.connect(self.open_semiexp_report)
        buttons.addWidget(open_report_button)
        self.semiexp_run_button = QPushButton("Run Semiexperimental Fit")
        self.semiexp_run_button.clicked.connect(self.run_semiexp_fit)
        buttons.addWidget(self.semiexp_run_button)
        run_layout.addLayout(buttons)
        tabs.addTab(input_tab, "Input")
        tabs.addTab(preview_tab, "GIC Preview")
        tabs.addTab(option_tab, "Fit Options")
        tabs.addTab(run_tab, "Run/Reports")
        self.add_semiexp_isotopologue_row(label="parent")
        return panel

    def semiexp_command_args(self) -> list[str]:
        args = [
            "semiexp",
            "--xyz",
            self.semiexp_xyz.text().strip(),
            "--observations",
            self.semiexp_observations.text().strip(),
            "--outdir",
            self.semiexp_outdir.text().strip(),
            "--backend",
            self.selected_backends.get("semiexp_geometry", "python"),
            "--observable",
            self.semiexp_observable.currentText(),
            "--rotational-components",
            self.semiexp_components.currentText(),
            "--prune-condition",
            f"{self.semiexp_prune_condition.value():.12g}",
        ]
        if self.semiexp_fixed.text().strip():
            args.extend(["--fixed", self.semiexp_fixed.text().strip()])
        for predicate in _split_semiexp_items(self.semiexp_qm.text()):
            args.extend(["--qm-predicate", predicate])
        for parameter_class in _split_semiexp_items(self.semiexp_classes.text()):
            args.extend(["--parameter-class", parameter_class])
        return args

    def _update_semiexp_preview(self) -> None:
        if not hasattr(self, "semiexp_command"):
            return
        args = self.semiexp_command_args()
        complete = all(args[idx] for idx in (2, 4, 6))
        command = "python -m merlino " + " ".join(shlex.quote(item) for item in args)
        if not complete:
            command += "\n\nSelect parent XYZ, isotopologue observations and output directory before running."
        self.semiexp_command.setPlainText(command)
        self.semiexp_run_button.setEnabled(complete)

    def run_semiexp_fit(self) -> None:
        args = self.semiexp_command_args()
        if not all(args[idx] for idx in (2, 4, 6)):
            self._update_semiexp_preview()
            return
        self._semiexp_process = QProcess(self)
        self._semiexp_process.setWorkingDirectory(str(self.workdir))
        self._semiexp_process.setProgram("python")
        self._semiexp_process.setArguments(["-m", "merlino", *args])
        self._semiexp_process.readyReadStandardOutput.connect(self._append_semiexp_output)
        self._semiexp_process.readyReadStandardError.connect(self._append_semiexp_output)
        self._semiexp_process.finished.connect(lambda code, status: self._append_semiexp_text(f"\nfinished: {code}\n"))
        self._append_semiexp_text("\nrunning...\n")
        self._semiexp_process.start()

    def add_semiexp_isotopologue_row(self, checked: bool = False, *, label: str = "") -> None:
        row = self.semiexp_iso_table.rowCount()
        self.semiexp_iso_table.insertRow(row)
        defaults = [label or f"iso{row + 1}", "", "", "", "", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "", "", ""]
        for col, value in enumerate(defaults):
            self.semiexp_iso_table.setItem(row, col, QTableWidgetItem(value))

    def save_semiexp_observations_toml(self, checked: bool = False) -> Path:
        target_text = self.semiexp_observations.text().strip()
        target = Path(target_text) if target_text else self.workdir / "isotopologues.toml"
        if target.suffix.lower() != ".toml":
            target = target.with_suffix(".toml")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(self._semiexp_table_toml(), encoding="utf-8")
        self.semiexp_observations.setText(str(target))
        self._append_semiexp_text(f"\nwrote observations: {target}\n")
        return target

    def preview_semiexp_gics(self) -> None:
        xyz = self.semiexp_xyz.text().strip()
        if not xyz:
            self._append_semiexp_text("\nselect parent XYZ before GIC preview\n")
            return
        observations = ()
        obs_path = self.semiexp_observations.text().strip()
        if obs_path and Path(obs_path).exists():
            observations = read_observations(Path(obs_path))
        preview = preview_semiexperimental_gics(Path(xyz), observations)
        self._fill_semiexp_preview_table(preview)
        self._append_semiexp_text("\n" + preview.text + "\n")

    def preview_semiexp_conditioning(self) -> None:
        request = self._semiexp_request_from_ui()
        if request is None:
            return
        preview = preview_semiexperimental_conditioning(request)
        self._append_semiexp_text("\n" + preview.text + "\n")

    def validate_semiexp_input(self) -> None:
        request = self._semiexp_request_from_ui()
        if request is None:
            return
        issues = validate_semiexperimental_request(request)
        if not issues:
            self._append_semiexp_text("\ninput validation: OK\n")
            return
        self._append_semiexp_text("\ninput validation:\n" + "\n".join(f"  {item.severity}: {item.message}" for item in issues) + "\n")

    def suggest_semiexp_classes(self) -> None:
        xyz = self.semiexp_xyz.text().strip()
        if not xyz:
            self._append_semiexp_text("\nselect parent XYZ before suggesting classes\n")
            return
        observations = ()
        obs_path = self.semiexp_observations.text().strip()
        if obs_path and Path(obs_path).exists():
            observations = read_observations(Path(obs_path))
        preview = preview_semiexperimental_gics(Path(xyz), observations)
        self._fill_semiexp_preview_table(preview)
        text = ";".join(f"{item.name}:{item.mode}:{'|'.join(item.patterns)}" for item in preview.suggested_classes)
        self.semiexp_classes.setText(text)
        self._append_semiexp_text("\n" + preview.text + "\n")

    def save_semiexp_preset(self, checked: bool = False) -> Path:
        path = self.workdir / "semiexp_preset.json"
        data = {
            "xyz": self.semiexp_xyz.text().strip(),
            "observations": self.semiexp_observations.text().strip(),
            "outdir": self.semiexp_outdir.text().strip(),
            "observable": self.semiexp_observable.currentText(),
            "components": self.semiexp_components.currentText(),
            "fixed": self.semiexp_fixed.text().strip(),
            "qm": self.semiexp_qm.text().strip(),
            "classes": self.semiexp_classes.text().strip(),
            "backend": self.selected_backends.get("semiexp_geometry", "python"),
        }
        path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        self._append_semiexp_text(f"\nwrote preset: {path}\n")
        return path

    def load_semiexp_preset(self, path: Path | bool | None = None) -> None:
        target = path if isinstance(path, Path) else self.workdir / "semiexp_preset.json"
        if not target.exists():
            self._append_semiexp_text(f"\npreset not found: {target}\n")
            return
        data = json.loads(target.read_text(encoding="utf-8"))
        self.semiexp_xyz.setText(str(data.get("xyz", "")))
        self.semiexp_observations.setText(str(data.get("observations", "")))
        self.semiexp_outdir.setText(str(data.get("outdir", "")))
        self.semiexp_observable.setCurrentText(str(data.get("observable", "moments")))
        self.semiexp_components.setCurrentText(str(data.get("components", "auto")))
        self.semiexp_fixed.setText(str(data.get("fixed", "")))
        self.semiexp_qm.setText(str(data.get("qm", "")))
        self.semiexp_classes.setText(str(data.get("classes", "")))
        backend = str(data.get("backend", "python"))
        if self.backend_selector.findText(backend) >= 0:
            self.backend_selector.setCurrentText(backend)
        self._append_semiexp_text(f"\nloaded preset: {target}\n")

    def open_semiexp_report(self) -> None:
        outdir = self.semiexp_outdir.text().strip()
        if not outdir:
            self._append_semiexp_text("\nselect output directory before opening report\n")
            return
        report = Path(outdir) / "semiexp_report.html"
        if report.exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(report)))
        else:
            self._append_semiexp_text(f"\nreport not found: {report}\n")

    def _append_semiexp_output(self) -> None:
        if self._semiexp_process is None:
            return
        text = bytes(self._semiexp_process.readAllStandardOutput()).decode(errors="replace")
        text += bytes(self._semiexp_process.readAllStandardError()).decode(errors="replace")
        self._append_semiexp_text(text)

    def _append_semiexp_text(self, text: str) -> None:
        self.semiexp_command.moveCursor(QTextCursor.MoveOperation.End)
        self.semiexp_command.insertPlainText(text)

    def _semiexp_table_toml(self) -> str:
        lines: list[str] = []
        for row in range(self.semiexp_iso_table.rowCount()):
            values = [_table_text(self.semiexp_iso_table, row, col) for col in range(self.semiexp_iso_table.columnCount())]
            if not values[0].strip():
                continue
            lines.extend([
                "[[isotopologues]]",
                f'label = "{_toml_string(values[0])}"',
                f'substitutions = "{_toml_string(values[1])}"',
                "[isotopologues.constants]",
                f"A_MHz = {_float_text(values[2])}",
                f"B_MHz = {_float_text(values[3])}",
                f"C_MHz = {_float_text(values[4])}",
                "[isotopologues.vibrational_correction]",
                f"delta_A_MHz = {_float_text(values[5], default='0.0')}",
                f"delta_B_MHz = {_float_text(values[6], default='0.0')}",
                f"delta_C_MHz = {_float_text(values[7], default='0.0')}",
                'source = "gui"',
                "[isotopologues.electronic_correction]",
                f"delta_A_MHz = {_float_text(values[8], default='0.0')}",
                f"delta_B_MHz = {_float_text(values[9], default='0.0')}",
                f"delta_C_MHz = {_float_text(values[10], default='0.0')}",
                'source = "gui"',
            ])
            if values[11].strip() and values[12].strip() and values[13].strip():
                lines.extend([
                    "[isotopologues.sigma_MHz]",
                    f"A_MHz = {_float_text(values[11])}",
                    f"B_MHz = {_float_text(values[12])}",
                    f"C_MHz = {_float_text(values[13])}",
                ])
            lines.append("")
        return "\n".join(lines)

    def _fill_semiexp_preview_table(self, preview) -> None:
        self.semiexp_preview_table.setRowCount(0)
        for row_data in preview.rows:
            row = self.semiexp_preview_table.rowCount()
            self.semiexp_preview_table.insertRow(row)
            values = [
                row_data.label,
                row_data.kind,
                ",".join(str(item) for item in row_data.atoms),
                row_data.suggested_class,
                row_data.state,
            ]
            for col, value in enumerate(values):
                self.semiexp_preview_table.setItem(row, col, QTableWidgetItem(value))

    def _semiexp_request_from_ui(self) -> SemiexperimentalFitRequest | None:
        xyz = self.semiexp_xyz.text().strip()
        obs_path = self.semiexp_observations.text().strip()
        if not xyz or not obs_path:
            self._append_semiexp_text("\nselect parent XYZ and observations before this operation\n")
            return None
        if not Path(obs_path).exists():
            self._append_semiexp_text(f"\nobservations file not found: {obs_path}\n")
            return None
        classes = []
        for item in _split_semiexp_items(self.semiexp_classes.text()):
            parts = item.split(":", 2)
            if len(parts) == 3:
                classes.append(ParameterClassConstraint(parts[0], tuple(p for p in parts[2].split("|") if p), parts[1]))
        return SemiexperimentalFitRequest(
            Path(xyz),
            read_observations(Path(obs_path)),
            fixed_parameters=tuple(_split_semiexp_items(self.semiexp_fixed.text().replace(",", ";"))),
            qm_predicates=tuple(_parse_semiexp_qm_predicates(self.semiexp_qm.text())),
            observable=self.semiexp_observable.currentText(),
            rotational_components=self.semiexp_components.currentText(),
            parameter_classes=tuple(classes),
        )


def workflow_detail_text(workflow: WorkflowSpec, selected_backend: str | None = None, workdir: Path | None = None) -> str:
    backend = selected_backend or workflow.default_backend
    lines = [
        workflow.title,
        "",
        workflow.description,
        "",
        "Category:",
        f"  {workflow.category}",
        "",
        "Service:",
        f"  {workflow.service}",
        "",
        "Backend:",
        f"  selected: {backend}",
        f"  available: {', '.join(workflow.backends)}",
        "",
        "Inputs:",
        *[f"  - {item}" for item in workflow.inputs],
        "",
        "Outputs:",
        *[f"  - {item}" for item in workflow.outputs],
        "",
        "Status:",
        f"  {workflow.status}",
    ]
    if workdir is not None:
        lines.extend(["", "Project State:", *[f"  - {item}" for item in workflow_state_lines(workflow, Path(workdir))]])
    return "\n".join(lines)


def workflow_state_lines(workflow: WorkflowSpec, workdir: Path) -> list[str]:
    manifests = _workflow_manifests(Path(workdir), workflow.workflow_id)
    if manifests:
        latest = manifests[-1]
        status = latest.get("status", "unknown")
        run_dir = latest.get("run_dir", "")
        lines = [f"latest manifest: {status} {run_dir}"]
        outputs = latest.get("outputs", {})
        if outputs:
            existing = [name for name, path in outputs.items() if Path(path).exists()]
            missing = [name for name, path in outputs.items() if not Path(path).exists()]
            lines.append(f"outputs present: {', '.join(existing) if existing else 'none'}")
            if missing:
                lines.append(f"outputs missing: {', '.join(missing)}")
        return lines
    expected = _expected_workflow_files(workflow.workflow_id, Path(workdir))
    existing = [name for name, path in expected.items() if path.exists()]
    if existing:
        return [f"files present without manifest: {', '.join(existing)}"]
    return ["not started in this workdir"]


def _workflow_manifests(workdir: Path, workflow_id: str) -> list[dict]:
    aliases = {
        "semiexp_geometry": {"semiexperimental_geometry"},
        "gic": {"gicforge"},
        "vpt2_vci": {"vpt2_vci", "gf"},
        "dvr": {"dvr"},
    }
    accepted = {workflow_id, *aliases.get(workflow_id, set())}
    data = []
    for path in sorted(workdir.rglob("*manifest*.json")):
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if manifest.get("workflow") in accepted:
            data.append(manifest)
    return data


def _expected_workflow_files(workflow_id: str, workdir: Path) -> dict[str, Path]:
    if workflow_id == "semiexp_geometry":
        return {
            "geometry": workdir / "semiexp_geometry.xyz",
            "parameters": workdir / "semiexp_parameters.csv",
            "geometry_parameters": workdir / "semiexp_geometry_parameters.csv",
            "report": workdir / "semiexp_report.html",
        }
    if workflow_id == "gic":
        return {"gaussian_input": workdir / "gauin", "report": workdir / "gicforge.out"}
    return {}


def _group_workflows(workflows: list[WorkflowSpec]) -> dict[str, list[WorkflowSpec]]:
    grouped: dict[str, list[WorkflowSpec]] = {}
    for workflow in workflows:
        grouped.setdefault(workflow.category, []).append(workflow)
    return grouped


def _path_row(parent: QWidget, form: QFormLayout, label: str, title: str, name_filter: str) -> QLineEdit:
    row = QHBoxLayout()
    edit = QLineEdit()
    edit.textChanged.connect(lambda _text: parent._update_semiexp_preview())
    button = QPushButton("Browse")
    button.clicked.connect(lambda: _select_file(parent, edit, title, name_filter))
    row.addWidget(edit, stretch=1)
    row.addWidget(button)
    form.addRow(label, row)
    return edit


def _directory_row(parent: QWidget, form: QFormLayout, label: str) -> QLineEdit:
    row = QHBoxLayout()
    edit = QLineEdit()
    edit.textChanged.connect(lambda _text: parent._update_semiexp_preview())
    button = QPushButton("Browse")
    button.clicked.connect(lambda: _select_directory(parent, edit))
    row.addWidget(edit, stretch=1)
    row.addWidget(button)
    form.addRow(label, row)
    return edit


def _select_file(parent: QWidget, edit: QLineEdit, title: str, name_filter: str) -> None:
    path, _selected_filter = QFileDialog.getOpenFileName(parent, title, str(parent.workdir), name_filter)
    if path:
        edit.setText(path)


def _select_directory(parent: QWidget, edit: QLineEdit) -> None:
    path = QFileDialog.getExistingDirectory(parent, "Select output directory", str(parent.workdir))
    if path:
        edit.setText(path)


def _split_semiexp_items(text: str) -> list[str]:
    return [item.strip() for item in text.split(";") if item.strip()]


def _parse_semiexp_qm_predicates(text: str) -> list[QMParameterPredicate]:
    predicates = []
    for item in _split_semiexp_items(text):
        parts = item.split(":")
        if len(parts) < 3:
            continue
        source = parts[3] if len(parts) > 3 else "gui"
        try:
            predicates.append(QMParameterPredicate(parts[0], float(parts[1]), float(parts[2]), source=source))
        except ValueError:
            continue
    return predicates


def _table_text(table: QTableWidget, row: int, col: int) -> str:
    item = table.item(row, col)
    return item.text().strip() if item is not None else ""


def _toml_string(text: str) -> str:
    return text.replace("\\", "\\\\").replace('"', '\\"')


def _float_text(text: str, *, default: str | None = None) -> str:
    raw = text.strip()
    if not raw and default is not None:
        return default
    float(raw)
    return raw
