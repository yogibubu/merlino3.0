from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QCheckBox,
    QMessageBox,
)
from PySide6.QtCore import Qt

from .gaussian import read_gaussian_alpha_with_freq
from .xyzin_utils import set_dvib_in_rotational


class DeltaVibAlphaDialog(QDialog):
    def __init__(self, log_path: Path, parent=None):
        super().__init__(parent)
        self.setWindowTitle("DeltaVib from Vibro-Rot alpha")
        self.resize(760, 520)

        self._log_path = Path(log_path)
        self._rows = []
        self._freq = []

        layout = QVBoxLayout(self)

        info = QLabel(
            "Compute ΔVib from Vibro-Rot alpha matrix (MHz). "
            "Default: invert alpha for imaginary frequencies."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        self.invert_imag_check = QCheckBox(
            "Invert alpha sign for imaginary frequencies"
        )
        self.invert_imag_check.setChecked(True)
        self.invert_imag_check.stateChanged.connect(self._refresh_table)
        layout.addWidget(self.invert_imag_check)

        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["Use", "Mode", "Freq (cm^-1)", "Imag", "Alpha A", "Alpha B", "Alpha C"]
        )
        self.table.itemChanged.connect(self._on_item_changed)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(self.table.SelectRows)
        layout.addWidget(self.table)

        sum_row = QHBoxLayout()
        self.sum_label = QLabel("Sum alpha (MHz): A=0.00000  B=0.00000  C=0.00000")
        self.dvib_label = QLabel("ΔVib (MHz): A=0.00000  B=0.00000  C=0.00000")
        self._dvib = (0.0, 0.0, 0.0)
        sum_row.addWidget(self.sum_label)
        sum_row.addStretch()
        sum_row.addWidget(self.dvib_label)
        layout.addLayout(sum_row)

        controls = QHBoxLayout()
        self.select_all_btn = QPushButton("Select all")
        self.select_none_btn = QPushButton("Select none")
        self.apply_btn = QPushButton("Apply to xyzin")
        self.close_btn = QPushButton("Close")

        self.select_all_btn.clicked.connect(lambda: self._set_all_checks(True))
        self.select_none_btn.clicked.connect(lambda: self._set_all_checks(False))
        self.apply_btn.clicked.connect(self._apply_to_xyzin)
        self.close_btn.clicked.connect(self.accept)

        controls.addWidget(self.select_all_btn)
        controls.addWidget(self.select_none_btn)
        controls.addStretch()
        controls.addWidget(self.apply_btn)
        controls.addWidget(self.close_btn)
        layout.addLayout(controls)

        self._load_data()

    def _load_data(self):
        rows, freq = read_gaussian_alpha_with_freq(self._log_path)
        if not rows:
            QMessageBox.warning(
                self,
                "DeltaVib",
                "Vibro-Rot alpha Matrix not found in Gaussian log.",
            )
            return

        self._rows = rows
        self._freq = freq or []
        self._refresh_table()

    def _freq_for_mode(self, mode: int):
        if mode - 1 < len(self._freq):
            return self._freq[mode - 1]
        return None

    def _effective_alpha(self, mode: int, a: float, b: float, c: float):
        f = self._freq_for_mode(mode)
        if self.invert_imag_check.isChecked() and f is not None and f < 0:
            return -a, -b, -c
        return a, b, c

    def _refresh_table(self):
        checked_modes = set()
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            mode_item = self.table.item(r, 1)
            if item is None or mode_item is None:
                continue
            if item.checkState() == Qt.Checked:
                try:
                    checked_modes.add(int(mode_item.text()))
                except Exception:
                    pass

        self.table.blockSignals(True)
        self.table.setRowCount(0)
        for row in self._rows:
            mode = row["mode"]
            f = self._freq_for_mode(mode)
            a_eff, b_eff, c_eff = self._effective_alpha(
                mode, row["a"], row["b"], row["c"]
            )
            r = self.table.rowCount()
            self.table.insertRow(r)

            item_use = QTableWidgetItem()
            item_use.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            if checked_modes and mode not in checked_modes:
                item_use.setCheckState(Qt.Unchecked)
            else:
                item_use.setCheckState(Qt.Checked)
            self.table.setItem(r, 0, item_use)

            self.table.setItem(r, 1, QTableWidgetItem(str(mode)))
            freq_txt = "" if f is None else f"{f:.6f}"
            self.table.setItem(r, 2, QTableWidgetItem(freq_txt))
            imag_txt = "Yes" if (f is not None and f < 0) else "No"
            self.table.setItem(r, 3, QTableWidgetItem(imag_txt))
            self.table.setItem(r, 4, QTableWidgetItem(f"{a_eff:.6f}"))
            self.table.setItem(r, 5, QTableWidgetItem(f"{b_eff:.6f}"))
            self.table.setItem(r, 6, QTableWidgetItem(f"{c_eff:.6f}"))

            for c in range(1, 7):
                item = self.table.item(r, c)
                if item is not None:
                    item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)

        self.table.blockSignals(False)
        self._recompute()

    def _set_all_checks(self, checked: bool):
        self.table.blockSignals(True)
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None:
                item.setCheckState(Qt.Checked if checked else Qt.Unchecked)
        self.table.blockSignals(False)
        self._recompute()

    def _on_item_changed(self, item: QTableWidgetItem):
        if item.column() == 0:
            self._recompute()

    def _recompute(self):
        sum_a = sum_b = sum_c = 0.0
        for r, row in enumerate(self._rows):
            item = self.table.item(r, 0)
            if item is None or item.checkState() != Qt.Checked:
                continue
            mode = row["mode"]
            a_eff, b_eff, c_eff = self._effective_alpha(
                mode, row["a"], row["b"], row["c"]
            )
            sum_a += a_eff
            sum_b += b_eff
            sum_c += c_eff

        self.sum_label.setText(
            f"Sum alpha (MHz): A={sum_a:.6f}  B={sum_b:.6f}  C={sum_c:.6f}"
        )
        self._dvib = (sum_a / 2, sum_b / 2, sum_c / 2)
        self.dvib_label.setText(
            f"ΔVib (MHz): A={self._dvib[0]:.6f}  B={self._dvib[1]:.6f}  C={self._dvib[2]:.6f}"
        )

    def _apply_to_xyzin(self):
        a, b, c = self._dvib
        ok = set_dvib_in_rotational(a, b, c)
        if not ok:
            QMessageBox.warning(self, "DeltaVib", "xyzin not found.")
            return
        QMessageBox.information(self, "DeltaVib", "Updated ΔVib in #ROTATIONAL.")
