from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMenu,
    QPushButton,
    QSplitter,
    QTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
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
        self.detail_view.setPlainText(workflow_detail_text(workflow, selected_backend=backend))

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
        self.detail_view.setPlainText(workflow_detail_text(workflow, selected_backend=backend))


def workflow_detail_text(workflow: WorkflowSpec, selected_backend: str | None = None) -> str:
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
    return "\n".join(lines)


def _group_workflows(workflows: list[WorkflowSpec]) -> dict[str, list[WorkflowSpec]]:
    grouped: dict[str, list[WorkflowSpec]] = {}
    for workflow in workflows:
        grouped.setdefault(workflow.category, []).append(workflow)
    return grouped
