from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QListWidget, QListWidgetItem, QMainWindow, QSplitter, QTextEdit, QVBoxLayout, QWidget


@dataclass(frozen=True)
class ManifestEntry:
    path: Path
    workflow: str
    status: str
    schema_version: str


def discover_manifests(root: Path) -> list[ManifestEntry]:
    """Discover Merlino manifest JSON files below a work directory."""
    entries: list[ManifestEntry] = []
    for path in sorted(Path(root).rglob("*manifest*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        entries.append(
            ManifestEntry(
                path=path,
                workflow=str(data.get("workflow", data.get("legacy", {}).get("workflow", "unknown"))),
                status=str(data.get("status", "unknown")),
                schema_version=str(data.get("schema_version", "legacy")),
            )
        )
    return entries


class ManifestBrowserWindow(QMainWindow):
    """Small browser for workflow manifests and reproducibility metadata."""

    def __init__(self, workdir: Path, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.entries = discover_manifests(self.workdir)
        self.setWindowTitle("Merlino Manifest Browser")
        self.resize(900, 600)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        layout = QVBoxLayout(central)
        self.setCentralWidget(central)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        self.list_widget = QListWidget()
        splitter.addWidget(self.list_widget)
        self.detail = QTextEdit()
        self.detail.setReadOnly(True)
        splitter.addWidget(self.detail)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)

        for entry in self.entries:
            item = QListWidgetItem(f"{entry.workflow} [{entry.status}] {entry.path.name}")
            item.setData(Qt.UserRole, str(entry.path))
            self.list_widget.addItem(item)
        self.list_widget.currentItemChanged.connect(self._show_manifest)
        if self.list_widget.count():
            self.list_widget.setCurrentRow(0)

    def _show_manifest(self, current: QListWidgetItem | None, previous=None) -> None:
        if current is None:
            self.detail.clear()
            return
        path = Path(current.data(Qt.UserRole))
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            self.detail.setPlainText(json.dumps(data, indent=2, sort_keys=True))
        except Exception as exc:
            self.detail.setPlainText(f"Could not read {path}: {exc}")
