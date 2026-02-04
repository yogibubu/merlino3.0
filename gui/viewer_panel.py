from pathlib import Path
import shutil
import subprocess
import tempfile
import os

from PySide6.QtWidgets import (
    QWidget, QLabel, QVBoxLayout, QPushButton,
    QHBoxLayout, QMessageBox, QMenuBar, QMenu
)
from PySide6.QtGui import QPixmap, QImage
from PySide6.QtCore import Qt

from . import viewer2d
from .xyz_reader import read_xyz


class ViewerPanel(QWidget):
    """
    Viewer panel for Merlino 3.0.

    - Shows Merlino logo by default
    - Renders molecule via viewer2d if xyzin is present
    - Avogadro (view) uses persistent temp file
    - Avogadro (edit) uses controlled temporary workspace
    """

    def __init__(self, working_dir: Path, parent=None):
        super().__init__(parent)

        self.working_dir = working_dir
        self.xyzin = self.working_dir / "xyzin"
        self._last_xyzin_mtime = None
        self._last_render_ok = False

        layout = QVBoxLayout(self)

        # --------------------------------------------------
        # Menu bar
        # --------------------------------------------------
        menubar = QMenuBar()
        view_menu = QMenu("View", self)
        menubar.addMenu(view_menu)

        action_view = view_menu.addAction("Open in Avogadro (view)")
        action_edit = view_menu.addAction("Open in Avogadro (edit)")
        action_edit_keep = view_menu.addAction(
            "Open in Avogadro (edit, keep SMILES)"
        )

        view_menu.addSeparator()

        action_view2 = view_menu.addAction("Open in Avogadro2 (view)")
        action_edit2 = view_menu.addAction("Open in Avogadro2 (edit)")

        action_view.triggered.connect(
            lambda: self._open_avogadro("avogadro", False, False)
        )
        action_edit.triggered.connect(
            lambda: self._open_avogadro("avogadro", True, False)
        )
        action_edit_keep.triggered.connect(
            lambda: self._open_avogadro("avogadro", True, True)
        )

        action_view2.triggered.connect(
            lambda: self._open_avogadro("avogadro2", False, False)
        )
        action_edit2.triggered.connect(
            lambda: self._open_avogadro("avogadro2", True, False)
        )

        layout.setMenuBar(menubar)

        # --------------------------------------------------
        # Viewer / logo
        # --------------------------------------------------
        self.label = QLabel(alignment=Qt.AlignCenter)
        layout.addWidget(self.label)

        # >>> DEFAULT LOGO (startup: 26)
        self.logo_path = Path(__file__).parent / "logo_prep.png"

        # --------------------------------------------------
        # Buttons
        # --------------------------------------------------
        btn_layout = QHBoxLayout()
        layout.addLayout(btn_layout)

        view_btn = QPushButton("Avogadro (view)")
        edit_btn = QPushButton("Avogadro (edit)")
        keep_btn = QPushButton("Edit (keep SMILES)")

        btn_layout.addWidget(view_btn)
        btn_layout.addWidget(edit_btn)
        btn_layout.addWidget(keep_btn)

        view_btn.clicked.connect(
            lambda: self._open_avogadro("avogadro", False, False)
        )
        edit_btn.clicked.connect(
            lambda: self._open_avogadro("avogadro", True, False)
        )
        keep_btn.clicked.connect(
            lambda: self._open_avogadro("avogadro", True, True)
        )

        self.show_logo()

    # --------------------------------------------------
    # Logo handling (NEW, MINIMAL)
    # --------------------------------------------------
    def set_logo(self, path: Path):
        """Change the logo shown when no molecule is displayed."""
        if path.exists():
            self.logo_path = path

    def show_logo(self):
        if self.logo_path.exists():
            pixmap = QPixmap(str(self.logo_path))
            self.label.setPixmap(
                pixmap.scaledToWidth(260, Qt.SmoothTransformation)
            )
            self.label.setStyleSheet(
                "background: #f7f3ee; border: 1px solid #eee6dd; "
                "padding: 24px 12px 12px 12px;"
            )
            self.label.setAlignment(Qt.AlignHCenter | Qt.AlignTop)
        else:
            self.label.setText("Merlino")

    # --------------------------------------------------
    # Viewer logic
    # --------------------------------------------------
    def refresh(self):
        if not self.xyzin.exists():
            self.show_logo()
            self._last_xyzin_mtime = None
            self._last_render_ok = False
            return

        try:
            mtime = self.xyzin.stat().st_mtime
        except OSError:
            self.show_logo()
            self._last_xyzin_mtime = None
            self._last_render_ok = False
            return

        if self._last_xyzin_mtime is not None and mtime == self._last_xyzin_mtime:
            if self._last_render_ok:
                return
        self._last_xyzin_mtime = mtime

        img = viewer2d.image_from_xyzin(self.xyzin)
        if img is None:
            self.show_logo()
            self._last_render_ok = False
            return

        img = img.convert("RGB")
        qimg = QImage(
            img.tobytes("raw", "RGB"),
            img.width,
            img.height,
            QImage.Format_RGB888,
        )

        pixmap = QPixmap.fromImage(qimg)
        self.label.setPixmap(
            pixmap.scaledToWidth(300, Qt.SmoothTransformation)
        )
        self.label.setStyleSheet("")
        self.label.setAlignment(Qt.AlignCenter)
        self._last_render_ok = True

    # --------------------------------------------------
    # Avogadro integration
    # --------------------------------------------------
    def _open_avogadro(self, cmd, import_back, keep_smiles):
        if not self.xyzin.exists():
            return

        if not import_back:
            fd, tmp_path = tempfile.mkstemp(suffix=".xyz")
            os.close(fd)
            shutil.copy(self.xyzin, tmp_path)

            try:
                if cmd == "avogadro":
                    subprocess.Popen(["open", tmp_path])
                else:
                    subprocess.Popen([cmd, tmp_path])
            except FileNotFoundError:
                QMessageBox.critical(
                    self,
                    "Executable not found",
                    f"'{cmd}' not found in PATH.",
                )
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_xyz = Path(tmpdir) / "structure.xyz"
            shutil.copy(self.xyzin, tmp_xyz)

            try:
                if cmd == "avogadro":
                    proc = subprocess.Popen(["open", str(tmp_xyz)])
                else:
                    proc = subprocess.Popen([cmd, str(tmp_xyz)])
            except FileNotFoundError:
                QMessageBox.critical(
                    self,
                    "Executable not found",
                    f"'{cmd}' not found in PATH.",
                )
                return

            proc.wait()

            reply = QMessageBox.question(
                self,
                "Import geometry",
                "Import modified geometry back into Merlino?",
                QMessageBox.Yes | QMessageBox.No,
            )

            if reply == QMessageBox.Yes:
                read_xyz(tmp_xyz)
                if keep_smiles:
                    self._restore_smiles()
                self.refresh()

    # --------------------------------------------------
    # Helpers
    # --------------------------------------------------
    def _restore_smiles(self):
        lines = self.xyzin.read_text().splitlines()
        smiles = []
        in_smiles = False
        for line in lines:
            if line.strip().upper() == "#SMILES":
                in_smiles = True
                continue
            if in_smiles:
                if line.startswith("#"):
                    break
                smiles.append(line.strip())

        if smiles:
            with self.xyzin.open("a") as f:
                f.write("#SMILES\n")
                for s in smiles:
                    f.write(s + "\n")
