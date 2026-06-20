from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import (
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from merlino_vpt2_vci import (
    VCIOptions,
    anharmonic_input_from_gaussian_fchk,
    compare_vpt2_vci,
    force_field_from_anharmonic_input,
    gf_from_hessian_input_with_merlino_gics,
    hessian_input_from_gaussian_fchk,
    read_indexed_qff_text,
)


class VPT2VCIWindow(QMainWindow):
    """Dedicated GUI for GF/PED and anharmonic VPT2/VCI workflows."""

    def __init__(self, workdir: Path, repo_root: Path | None = None, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.repo_root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]

        self.setWindowTitle("Merlino GF / VPT2-VCI")
        self.resize(980, 780)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        header = QLabel("GF / VPT2-VCI")
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        note = QLabel(
            "Gaussian/FCHK is only an input adapter here. GF, PED, VPT2 and VCI "
            "run on canonical Merlino data structures and Merlino non-redundant GICs."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        input_group = QGroupBox("Inputs")
        input_layout = QVBoxLayout(input_group)

        row_fchk = QHBoxLayout()
        row_fchk.addWidget(QLabel("FCHK Hessian/QFF:"))
        self.fchk_edit = QLineEdit(str(self.workdir / "gauin.fchk"))
        row_fchk.addWidget(self.fchk_edit)
        browse_fchk = QPushButton("Browse")
        browse_fchk.clicked.connect(self._browse_fchk)
        row_fchk.addWidget(browse_fchk)
        input_layout.addLayout(row_fchk)

        row_qff = QHBoxLayout()
        row_qff.addWidget(QLabel("Indexed QFF text:"))
        self.qff_edit = QLineEdit("")
        row_qff.addWidget(self.qff_edit)
        browse_qff = QPushButton("Browse")
        browse_qff.clicked.connect(self._browse_qff)
        row_qff.addWidget(browse_qff)
        input_layout.addLayout(row_qff)

        row_latest = QHBoxLayout()
        latest_fchk = QPushButton("Use Latest FCHK")
        latest_fchk.clicked.connect(lambda: self._select_latest_fchk())
        row_latest.addWidget(latest_fchk)
        row_latest.addStretch()
        input_layout.addLayout(row_latest)
        layout.addWidget(input_group)

        gf_group = QGroupBox("GF / PED")
        gf_layout = QVBoxLayout(gf_group)
        gf_note = QLabel(
            "Reads Cartesian Hessian and geometry from FCHK, builds Merlino GICs, "
            "then solves Wilson GF and reports PED in non-redundant GICs."
        )
        gf_note.setWordWrap(True)
        gf_layout.addWidget(gf_note)
        self.run_gf_button = QPushButton("Run GF / PED")
        self.run_gf_button.clicked.connect(self.run_gf)
        gf_layout.addWidget(self.run_gf_button)
        layout.addWidget(gf_group)

        vci_group = QGroupBox("VPT2 / VCI")
        vci_layout = QVBoxLayout(vci_group)

        row_basis = QHBoxLayout()
        row_basis.addWidget(QLabel("Max total quanta:"))
        self.max_quanta_edit = QLineEdit("2")
        row_basis.addWidget(self.max_quanta_edit)
        row_basis.addWidget(QLabel("Roots:"))
        self.roots_edit = QLineEdit("6")
        row_basis.addWidget(self.roots_edit)
        row_basis.addWidget(QLabel("Active modes (1-based):"))
        self.active_modes_edit = QLineEdit("")
        self.active_modes_edit.setPlaceholderText("blank = all, e.g. 1,2,5")
        row_basis.addWidget(self.active_modes_edit)
        vci_layout.addLayout(row_basis)

        row_prune = QHBoxLayout()
        row_prune.addWidget(QLabel("Freq min cm-1:"))
        self.freq_min_edit = QLineEdit("")
        row_prune.addWidget(self.freq_min_edit)
        row_prune.addWidget(QLabel("Freq max cm-1:"))
        self.freq_max_edit = QLineEdit("")
        row_prune.addWidget(self.freq_max_edit)
        row_prune.addWidget(QLabel("Basis cutoff cm-1:"))
        self.basis_cutoff_edit = QLineEdit("")
        row_prune.addWidget(self.basis_cutoff_edit)
        row_prune.addWidget(QLabel("Force threshold cm-1:"))
        self.force_threshold_edit = QLineEdit("0.0")
        row_prune.addWidget(self.force_threshold_edit)
        vci_layout.addLayout(row_prune)

        row_limits = QHBoxLayout()
        row_limits.addWidget(QLabel("Mode max quanta:"))
        self.mode_max_edit = QLineEdit("")
        self.mode_max_edit.setPlaceholderText("blank or per original mode, e.g. 3,3,2")
        row_limits.addWidget(self.mode_max_edit)
        row_limits.addWidget(QLabel("Class limits 1/2/3/4-mode:"))
        self.class_limits_edit = QLineEdit("")
        self.class_limits_edit.setPlaceholderText("e.g. 1:1-2;2:2-3;3:3-4;4:4-4")
        row_limits.addWidget(self.class_limits_edit)
        vci_layout.addLayout(row_limits)

        row_run = QHBoxLayout()
        self.run_vpt2_vci_button = QPushButton("Run VPT2 / VCI")
        self.run_vpt2_vci_button.clicked.connect(self.run_vpt2_vci)
        row_run.addWidget(self.run_vpt2_vci_button)
        row_run.addStretch()
        vci_layout.addLayout(row_run)
        layout.addWidget(vci_group)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setPlaceholderText("GF/PED and VPT2/VCI reports will appear here.")

        actions = QHBoxLayout()
        clear_button = QPushButton("Clear Output")
        clear_button.clicked.connect(self.output_text.clear)
        actions.addWidget(clear_button)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        actions.addWidget(close_button)
        actions.addStretch()
        layout.addLayout(actions)

        layout.addWidget(self.output_text)

    def run_gf(self, *, show_message: bool = True) -> None:
        try:
            fchk_path = self._required_existing_path(self.fchk_edit.text(), "FCHK")
            hessian_input = hessian_input_from_gaussian_fchk(fchk_path)
            result = gf_from_hessian_input_with_merlino_gics(hessian_input)
            self.output_text.setPlainText(self._format_gf_report(fchk_path, result))
        except Exception as exc:
            self._fail("GF / PED failed", exc, show_message)

    def run_vpt2_vci(self, *, show_message: bool = True) -> None:
        try:
            qff = self._load_force_field()
            max_quanta = self._parse_int(self.max_quanta_edit.text(), "Max total quanta", minimum=0)
            roots = self._parse_int(self.roots_edit.text(), "Roots", minimum=1)
            options = VCIOptions(
                active_modes=self._parse_active_modes(),
                frequency_min_cm=self._parse_optional_float(self.freq_min_edit.text(), "Freq min cm-1"),
                frequency_max_cm=self._parse_optional_float(self.freq_max_edit.text(), "Freq max cm-1"),
                basis_energy_cutoff_cm=self._parse_optional_float(self.basis_cutoff_edit.text(), "Basis cutoff cm-1"),
                mode_max_quanta=self._parse_optional_int_tuple(self.mode_max_edit.text(), "Mode max quanta"),
                excitation_class_limits=self._parse_class_limits(),
                force_constant_threshold_cm=self._parse_optional_float(
                    self.force_threshold_edit.text(), "Force threshold cm-1", default=0.0
                ),
            )
            comparison = compare_vpt2_vci(qff, max_quanta=max_quanta, n_roots=roots, options=options)
            self.output_text.setPlainText(self._format_vpt2_vci_report(qff, comparison))
        except Exception as exc:
            self._fail("VPT2 / VCI failed", exc, show_message)

    def _load_force_field(self):
        qff_path = self._optional_existing_path(self.qff_edit.text())
        fchk_path = self._optional_existing_path(self.fchk_edit.text())

        frequencies = None
        if fchk_path is not None:
            anharmonic_input = anharmonic_input_from_gaussian_fchk(fchk_path)
            frequencies = (
                anharmonic_input.anharmonic_frequencies_cm
                if anharmonic_input.anharmonic_frequencies_cm.size
                else anharmonic_input.harmonic_frequencies_cm
            )
        if qff_path is not None:
            return read_indexed_qff_text(qff_path, frequencies)
        if fchk_path is not None:
            return force_field_from_anharmonic_input(anharmonic_input)
        raise FileNotFoundError("Provide an existing FCHK file or indexed QFF text file")

    def _format_gf_report(self, fchk_path: Path, result) -> str:
        lines = [
            "GF/PED from Merlino non-redundant GICs",
            f"Source FCHK: {fchk_path}",
            f"GIC count: {len(result.gic_labels)}",
            "",
            "Frequencies (cm-1):",
        ]
        for idx, freq in enumerate(result.frequencies_cm, start=1):
            lines.append(f"  mode {idx:3d}: {freq:12.3f}")

        lines.extend(["", "GIC labels:"])
        for idx, label in enumerate(result.gic_labels, start=1):
            lines.append(f"  GIC{idx:03d}: {label}")

        lines.extend(["", "PED (%) rows=GIC cols=modes:"])
        header = "          " + " ".join(f"M{idx:02d}" for idx in range(1, len(result.frequencies_cm) + 1))
        lines.append(header)
        for idx, row in enumerate(result.ped.values, start=1):
            values = " ".join(f"{value:7.2f}" for value in row)
            lines.append(f"  GIC{idx:03d} {values}")
        return "\n".join(lines)

    def _format_vpt2_vci_report(self, qff, comparison) -> str:
        lines = [
            "VPT2/VCI comparison on canonical Merlino QFF",
            f"Modes used in input force field: {len(qff.harmonic_frequencies_cm)}",
            f"Cubic terms: {len(qff.cubic_cm)}",
            f"Quartic terms: {len(qff.quartic_cm)}",
            f"VCI basis size: {len(comparison.vci.basis)}",
            "Input harmonic frequencies (cm-1): "
            + ", ".join(f"{value:.3f}" for value in qff.harmonic_frequencies_cm),
            "",
            "Root     VPT2 abs      VCI abs        d_abs     VPT2 exc      VCI exc        d_exc",
        ]
        n = min(
            len(comparison.vpt2.energies_cm),
            len(comparison.vci.energies_cm),
            len(comparison.energy_differences_cm),
        )
        for idx in range(n):
            lines.append(
                f"{idx + 1:4d} "
                f"{comparison.vpt2.energies_cm[idx]:12.4f} "
                f"{comparison.vci.energies_cm[idx]:12.4f} "
                f"{comparison.energy_differences_cm[idx]:10.4f} "
                f"{comparison.vpt2.excitation_energies_cm[idx]:12.4f} "
                f"{comparison.vci.excitation_energies_cm[idx]:12.4f} "
                f"{comparison.excitation_differences_cm[idx]:10.4f}"
            )

        if comparison.vci.blocks:
            lines.extend(["", "Symmetry blocks:"])
            for block in comparison.vci.blocks:
                lines.append(f"  {block.label}: states={len(block.basis_indices)} roots={block.n_roots}")

        if comparison.vci.state_contributions:
            lines.extend(["", "Dominant VCI contributions:"])
            for root, contribution in enumerate(comparison.vci.state_contributions[:n], start=1):
                pieces = [
                    f"{state}:{coeff:+.3f}"
                    for state, coeff in contribution.dominant_basis_states[:4]
                ]
                lines.append(
                    f"  root {root:3d}: <n>={np.array2string(contribution.mode_quanta, precision=3)} "
                    + ", ".join(pieces)
                )
        return "\n".join(lines)

    def _browse_fchk(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select FCHK", str(self.workdir), "FCHK files (*.fchk *.fch);;All files (*)")
        if path:
            self.fchk_edit.setText(path)

    def _browse_qff(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select indexed QFF", str(self.workdir), "QFF text (*.qff *.txt);;All files (*)")
        if path:
            self.qff_edit.setText(path)

    def _select_latest_fchk(self, *, show_message: bool = True) -> None:
        candidates = sorted(
            list(self.workdir.glob("*.fchk")) + list(self.workdir.glob("*.fch")),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if not candidates:
            if show_message:
                QMessageBox.warning(self, "FCHK", f"No FCHK files found in {self.workdir}")
            return
        self.fchk_edit.setText(str(candidates[0]))

    def _required_existing_path(self, raw: str, label: str) -> Path:
        path = Path(raw.strip()).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
        return path

    def _optional_existing_path(self, raw: str) -> Path | None:
        text = raw.strip()
        if not text:
            return None
        path = Path(text).expanduser()
        if not path.exists():
            return None
        return path

    def _parse_active_modes(self) -> tuple[int, ...] | None:
        values = self._parse_optional_int_tuple(self.active_modes_edit.text(), "Active modes")
        if values is None:
            return None
        if any(value < 1 for value in values):
            raise ValueError("Active modes are one-based and must be positive")
        return tuple(value - 1 for value in values)

    def _parse_class_limits(self) -> dict[int, tuple[int, int | None]]:
        text = self.class_limits_edit.text().strip()
        if not text:
            return {}
        limits: dict[int, tuple[int, int | None]] = {}
        for item in text.replace(",", ";").split(";"):
            item = item.strip()
            if not item:
                continue
            if ":" not in item:
                raise ValueError("Class limits must use n:min-max records")
            key_raw, value_raw = item.split(":", 1)
            n_modes = self._parse_int(key_raw, "Class limit key", minimum=1)
            if "-" in value_raw:
                lo_raw, hi_raw = value_raw.split("-", 1)
                qmin = self._parse_int(lo_raw, "Class limit minimum", minimum=0)
                qmax = None if hi_raw.strip() == "*" else self._parse_int(hi_raw, "Class limit maximum", minimum=qmin)
            else:
                qmin = self._parse_int(value_raw, "Class limit minimum", minimum=0)
                qmax = qmin
            limits[n_modes] = (qmin, qmax)
        return limits

    @staticmethod
    def _parse_optional_int_tuple(raw: str, label: str) -> tuple[int, ...] | None:
        text = raw.strip()
        if not text:
            return None
        values = []
        for token in text.replace(";", ",").split(","):
            token = token.strip()
            if token:
                values.append(VPT2VCIWindow._parse_int(token, label, minimum=0))
        return tuple(values)

    @staticmethod
    def _parse_optional_float(raw: str, label: str, *, default: float | None = None) -> float | None:
        text = raw.strip()
        if not text:
            return default
        try:
            return float(text)
        except ValueError as exc:
            raise ValueError(f"{label} must be numeric") from exc

    @staticmethod
    def _parse_int(raw: str, label: str, *, minimum: int) -> int:
        try:
            value = int(raw.strip())
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer") from exc
        if value < minimum:
            raise ValueError(f"{label} must be >= {minimum}")
        return value

    def _fail(self, title: str, exc: Exception, show_message: bool) -> None:
        message = f"{title}: {exc}"
        self.output_text.setPlainText(message)
        if show_message:
            QMessageBox.critical(self, title, str(exc))
