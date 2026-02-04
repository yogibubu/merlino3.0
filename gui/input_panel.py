from pathlib import Path

from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QRadioButton,
    QButtonGroup,
    QMessageBox,
    QCheckBox,
)

from .readers import read_structure
from .gaussian import read_gaussian_properties


class InputSource:
    """Minimal DTO for readers."""

    def __init__(self, kind, path=None, smiles=None):
        self.kind = kind
        self.path = path
        self.smiles = smiles


class InputPanel(QWidget):
    def __init__(self, working_dir: Path, on_update, parent=None):
        super().__init__(parent)

        self.on_update = on_update
        self.working_dir = working_dir

        main_layout = QVBoxLayout(self)

        # ==================================================
        # Input type selector
        # ==================================================
        type_layout = QHBoxLayout()
        type_label = QLabel("Input type:")
        type_layout.addWidget(type_label)

        self._type_group = QButtonGroup(self)
        self._type_buttons = {}

        for i, name in enumerate(
            ["SMILES", "XYZ", "Z-matrix", "Gaussian", "Molpro", "MRCC"]
        ):
            btn = QRadioButton(name)
            self._type_group.addButton(btn, i)
            self._type_buttons[name] = btn
            type_layout.addWidget(btn)

        type_layout.addStretch()
        main_layout.addLayout(type_layout)

        self._type_buttons["SMILES"].setChecked(True)
        self._type_group.idClicked.connect(self._update_mode)

        # ==================================================
        # SMILES input
        # ==================================================
        smiles_layout = QHBoxLayout()
        smiles_label = QLabel("SMILES:")
        self.smiles_edit = QLineEdit()
        self.smiles_edit.setPlaceholderText(
            "Paste or type SMILES, press Enter to commit"
        )

        smiles_layout.addWidget(smiles_label)
        smiles_layout.addWidget(self.smiles_edit)
        main_layout.addLayout(smiles_layout)

        self.smiles_edit.returnPressed.connect(self._on_smiles_commit)

        # ==================================================
        # FILE input
        # ==================================================
        file_layout = QHBoxLayout()
        file_label = QLabel("File:")
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText(
            "Paste or type file path, press Enter to commit"
        )
        self.browse_btn = QPushButton("Browse…")

        file_layout.addWidget(file_label)
        file_layout.addWidget(self.file_edit)
        file_layout.addWidget(self.browse_btn)
        main_layout.addLayout(file_layout)

        self.file_edit.returnPressed.connect(self._on_file_commit)
        self.browse_btn.clicked.connect(self._on_browse_file)

        # ==================================================
        # Optional Gaussian properties file
        # ==================================================
        props_layout = QHBoxLayout()
        self.props_check = QCheckBox("Use Gaussian properties")
        self.props_check.setToolTip(
            "Gaussian .log/.out only. FCHK/FCH are used in other pipelines (normal modes / Hessian)."
        )
        self.props_edit = QLineEdit()
        self.props_edit.setPlaceholderText("Gaussian .log/.out for properties only")
        self.props_edit.setToolTip("Select a Gaussian .log/.out; geometry will NOT be replaced.")
        self.props_browse = QPushButton("Browse…")

        props_layout.addWidget(self.props_check)
        props_layout.addWidget(self.props_edit)
        props_layout.addWidget(self.props_browse)
        main_layout.addLayout(props_layout)

        self.props_check.stateChanged.connect(self._update_mode)
        self.props_edit.returnPressed.connect(self._on_file_commit)
        self.props_browse.clicked.connect(self._on_browse_props)

        # Initial enable/disable
        self._update_mode()

    # ==================================================
    # Utilities
    # ==================================================
    def _show_error(self, title: str, message: str):
        QMessageBox.critical(self, title, message)

    # ==================================================
    # Mode handling
    # ==================================================
    def _current_mode(self) -> str:
        btn = self._type_group.checkedButton()
        return btn.text() if btn else "SMILES"

    def _update_mode(self, *_):
        mode = self._current_mode()

        is_smiles = mode == "SMILES"
        is_file = mode in {"XYZ", "Z-matrix", "Gaussian", "Molpro", "MRCC"}

        self.smiles_edit.setEnabled(is_smiles)
        self.file_edit.setEnabled(is_file)
        self.browse_btn.setEnabled(is_file)

        if is_smiles:
            self.file_edit.clear()
        else:
            self.smiles_edit.clear()

        # Gaussian properties file can be provided for any mode
        allow_props = True
        self.props_check.setEnabled(allow_props)
        self.props_edit.setEnabled(allow_props and self.props_check.isChecked())
        self.props_browse.setEnabled(allow_props and self.props_check.isChecked())

    # ==================================================
    # SMILES handler (COMMIT ONLY)
    # ==================================================
    def _on_smiles_commit(self):
        if self._current_mode() != "SMILES":
            return

        text = self.smiles_edit.text().strip()
        if not text:
            return

        src = InputSource(kind="smiles", smiles=text)

        try:
            read_structure(src)
            if self.props_check.isChecked():
                if not self.props_edit.text().strip():
                    self._show_error("Input error", "Properties file not specified.")
                    return
                ppath = Path(self.props_edit.text().strip())
                if not ppath.exists():
                    self._show_error("Input error", f"Properties file not found:\n{ppath}")
                    return
                if ppath.suffix.lower() in {".fchk", ".fch"}:
                    self._show_error(
                        "Input error",
                        "FCHK/FCH are not used here. Use a Gaussian .log/.out for properties.",
                    )
                    return
                read_gaussian_properties(ppath)
        except Exception as e:
            self._show_error("Input error", str(e))
            return

        self.on_update("smiles")

    # ==================================================
    # FILE handlers
    # ==================================================
    def _on_file_commit(self):
        mode = self._current_mode()
        if mode == "SMILES":
            return

        path = Path(self.file_edit.text().strip())
        if not path.exists():
            self._show_error("Input error", f"File not found:\n{path}")
            return

        kind = "zmat" if mode == "Z-matrix" else mode.lower()
        src = InputSource(kind=kind, path=path)

        try:
            read_structure(src)
            if self.props_check.isChecked():
                if not self.props_edit.text().strip():
                    self._show_error("Input error", "Properties file not specified.")
                    return
                ppath = Path(self.props_edit.text().strip())
                if not ppath.exists():
                    self._show_error("Input error", f"Properties file not found:\n{ppath}")
                    return
                if ppath.suffix.lower() in {".fchk", ".fch"}:
                    self._show_error(
                        "Input error",
                        "FCHK/FCH are not used here. Use a Gaussian .log/.out for properties.",
                    )
                    return
                read_gaussian_properties(ppath)
        except Exception as e:
            self._show_error("Input error", str(e))
            return

        self.on_update(kind)

    def _on_browse_file(self):
        mode = self._current_mode()
        if mode not in {"XYZ", "Z-matrix", "Gaussian", "Molpro", "MRCC"}:
            return

        filter_map = {
            "XYZ": "XYZ files (*.xyz);;All files (*)",
            "Z-matrix": "Z-matrix files (*.zmat *.inp);;All files (*)",
            "Gaussian": "Gaussian output (*.log *.out *.fchk *.fch);;All files (*)",
            "Molpro": "Molpro output (*.out *.log);;All files (*)",
            "MRCC": "MRCC output (*.out *.log);;All files (*)",
        }

        file_filter = filter_map.get(mode, "All files (*)")

        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Select {mode} input file",
            "",
            file_filter,
        )
        if not path:
            return

        self.file_edit.setText(path)
        self._on_file_commit()

    def _on_browse_props(self):
        if not self.props_check.isChecked():
            return
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Gaussian properties file",
            "",
            "Gaussian output (*.log *.out);;All files (*)",
        )
        if not path:
            return
        self.props_edit.setText(path)
        self._on_file_commit()
