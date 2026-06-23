from __future__ import annotations

from pathlib import Path
import shutil

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QGroupBox,
    QMessageBox, QCheckBox, QComboBox,
    QScrollArea, QInputDialog, QLineEdit, QFileDialog
)
from PySide6.QtCore import Qt
import numpy as np
from PySide6.QtGui import QPixmap

from merlino_core import ensure_workspace, load_config

# Merlino / QM
from advanced.provin_writer import ProvinWriter
from advanced.launchers.gicforge_launcher import GICForgeLauncher
from advanced.launchers.msr_launcher import MSRLauncher
from advanced.launchers.gaussian_launcher import GaussianLauncher
from advanced.launchers.survibfit_launcher import SurvibfitLauncher
from advanced.dvr_window import DVRWindow
from advanced.gf_window import GFWindow
from advanced.vpt2_vci_window import VPT2VCIWindow

from advanced.kwd_spec import KWD_SPEC

# Viewers
from advanced.viewers.gicforge_viewer import GICForgeViewer
from advanced.viewers.gaussian_viewer import GaussianViewer


ADVANCED_SECTIONS = {
    "geometry_topology",
    "gnic",
    "method",
    "calculation",
    "symmetry",
    "writer",
    "workflow",
}


class AdvancedWindow(QMainWindow):
    """
    Merlino / QM – Advanced Calculations Window
    """

    def __init__(self, workdir: Path, parent=None):
        super().__init__(parent)

        self.workdir = Path(workdir)
        self.project_root = Path(__file__).resolve().parents[1]
        self.project_dir = self.workdir.parent
        self.xyzin = self.workdir / "xyzin"
        self.project_name = None
        self.project_out_dir = None
        self.config = load_config(workdir=self.workdir)
        self.workspace = None

        self.kwd_widgets = {}

        self.setWindowTitle("Merlino / QM – Advanced Calculations")
        self.resize(800, 600)

        self._init_project_dir()
        self._build_ui()

    # ==================================================================
    # UI
    # ==================================================================

    def _build_ui(self):
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        central = QWidget()
        scroll.setWidget(central)
        self.setCentralWidget(scroll)

        layout = QVBoxLayout(central)

        module_dir = Path(__file__).resolve().parent
        logo_path = module_dir / "logo_adv.png"

        if logo_path.exists():
            logo = QLabel()
            logo.setPixmap(
                QPixmap(str(logo_path)).scaledToWidth(
                    300, Qt.SmoothTransformation
                )
            )
            logo.setAlignment(Qt.AlignCenter)
            layout.addWidget(logo)

        header = QLabel("Advanced Calculations")
        header.setAlignment(Qt.AlignCenter)
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        kw_group = QGroupBox("Preparation keywords")
        kw_layout = QVBoxLayout(kw_group)

        for sec_name, section in KWD_SPEC.items():
            if sec_name not in ADVANCED_SECTIONS:
                continue

            box = QGroupBox(section["panel"])
            row = QHBoxLayout(box)

            if sec_name == "gnic":
                gnic_cb = QCheckBox("GNIC")
                gnic_cb.setChecked(True)
                gnic_cb.setEnabled(False)
                self.kwd_widgets["GNIC"] = gnic_cb
                row.addWidget(gnic_cb)

                for kw in ["NOONEDIH", "BDPCS3"]:
                    cb = QCheckBox(kw)
                    self.kwd_widgets[kw] = cb
                    row.addWidget(cb)

                kw_layout.addWidget(box)
                continue

            for kw, spec in section["keywords"].items():
                if spec.get("type") == "flag":
                    cb = QCheckBox(kw)
                    cb.setChecked(spec.get("default", False))
                    self.kwd_widgets[kw] = cb
                    row.addWidget(cb)
                elif "values" in spec:
                    combo = QComboBox()
                    for v in spec["values"]:
                        if v is not None:
                            combo.addItem(v)
                    if spec.get("default") is not None:
                        combo.setCurrentText(spec["default"])
                    self.kwd_widgets[kw] = combo
                    row.addWidget(QLabel(kw))
                    row.addWidget(combo)

            kw_layout.addWidget(box)

        layout.addWidget(kw_group)

        self.run_prova_btn = QPushButton("Run GICForge")
        self.run_prova_btn.clicked.connect(lambda: AdvancedWindow.run_gicforge(self))
        layout.addWidget(self.run_prova_btn)

        self.run_mfit_btn = QPushButton("Prepare gauin.gjf (merlino_fit)")
        self.run_mfit_btn.clicked.connect(lambda: AdvancedWindow.run_merlino_fit_prep(self))
        layout.addWidget(self.run_mfit_btn)

        methods_layout = QHBoxLayout()

        self.run_msr_btn = QPushButton("Run MSR")
        self.run_msr_btn.clicked.connect(lambda: AdvancedWindow.run_msr(self))
        self.run_msr_btn.setEnabled(False)

        self.run_gaussian_btn = QPushButton("Run Gaussian")
        self.run_gaussian_btn.clicked.connect(lambda: AdvancedWindow.run_gaussian(self))
        self.run_gaussian_btn.setEnabled(False)

        methods_layout.addWidget(self.run_msr_btn)
        methods_layout.addWidget(self.run_gaussian_btn)

        layout.addLayout(methods_layout)

        self._build_survibfit_panel(layout)
        self._build_gf_shortcut_panel(layout)
        self._build_vpt2_vci_shortcut_panel(layout)
        self._build_dvr_shortcut_panel(layout)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

    # ==================================================================
    # Helpers
    # ==================================================================
    def _init_project_dir(self):
        base = self.project_root / "projects"
        base.mkdir(parents=True, exist_ok=True)

        name, ok = QInputDialog.getText(
            self,
            "Project name",
            "Project name:",
        )
        name = (name or "").strip()
        if not ok or not name:
            name = self.workdir.name

        self.project_name = name
        self.project_out_dir = base / name
        self.project_out_dir.mkdir(parents=True, exist_ok=True)
        self.workspace = ensure_workspace(self.project_out_dir)

    def _get_merlino_fit_root(self) -> Path:
        root = self.project_root / "merlino_fit"
        if root.exists():
            return root
        return Path("/Users/vincenzobarone/merlino_fit")


    def _export_project_files(self):
        if not self.project_out_dir or not self.project_name:
            return

        mapping = {
            self.workdir / "xyzin": f"{self.project_name}.xyz",
            self.workdir / "gauin.gjf": f"{self.project_name}.gjf",
            self.workdir / "gauin": f"{self.project_name}.gjf",
            self.workdir / "gicforge.chk": f"{self.project_name}.chk",
            self.workdir / "gicforge.out": f"{self.project_name}.ncc",
        }

        log_candidates = []
        for log_name in ("gauin.log", "gauout.log"):
            p = self.workdir / log_name
            if p.exists():
                try:
                    st = p.stat()
                except Exception:
                    continue
                log_candidates.append((st.st_mtime, st.st_size, p))
        if log_candidates:
            log_candidates.sort(key=lambda t: (t[0], t[1]), reverse=True)
            mapping[log_candidates[0][2]] = f"{self.project_name}.log"

        for src, dst_name in mapping.items():
            if src.exists():
                dst = self.project_out_dir / dst_name
                try:
                    shutil.copy2(src, dst)
                except Exception:
                    continue

    def _read_charge_multiplicity(self):
        if not self.xyzin.exists():
            raise FileNotFoundError("xyzin not found")

        charge = None
        multiplicity = None
        in_basic = False

        for line in self.xyzin.read_text().splitlines():
            line = line.strip()
            if not line:
                continue

            if line.upper().startswith("#BASIC"):
                in_basic = True
                continue

            if in_basic and line.startswith("#"):
                break

            if in_basic:
                if "=" in line:
                    key, val = [x.strip() for x in line.split("=", 1)]
                    parts = [key, val]
                else:
                    parts = line.split()
                if len(parts) < 2:
                    continue
                key = parts[0].strip().lower()
                val = parts[-1].strip()
                if key == "charge":
                    charge = int(val)
                elif key in ("multiplicity", "spin_multiplicity"):
                    multiplicity = int(val)

        if charge is None or multiplicity is None:
            raise ValueError("Missing CHARGE or SPIN_MULTIPLICITY in #BASIC")

        return charge, multiplicity

    def _collect_gicforge_keywords(self) -> list[str]:
        keywords = []
        for kw, widget in self.kwd_widgets.items():
            if isinstance(widget, QCheckBox):
                if widget.isChecked():
                    keywords.append(kw)
            elif isinstance(widget, QComboBox):
                val = widget.currentText()
                if val:
                    keywords.append(val)
        return keywords

    def _collect_prova_keywords(self) -> list[str]:
        return self._collect_gicforge_keywords()

    # ==================================================================
    # Survibfit panel
    # ==================================================================
    def _build_survibfit_panel(self, layout: QVBoxLayout):
        group = QGroupBox("Survibfit – Vibrational analysis")
        vbox = QVBoxLayout(group)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Log file:"))
        self.sv_log = QLineEdit(str(self.workdir / "gauin.log"))
        row1.addWidget(self.sv_log)
        btn_log = QPushButton("Browse")
        btn_log.clicked.connect(self._browse_log)
        row1.addWidget(btn_log)
        vbox.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("FCHK file:"))
        self.sv_fchk = QLineEdit(str(self.workdir / "gauin.fchk"))
        row2.addWidget(self.sv_fchk)
        btn_fchk = QPushButton("Browse")
        btn_fchk.clicked.connect(self._browse_fchk)
        row2.addWidget(btn_fchk)
        vbox.addLayout(row2)

        row3 = QHBoxLayout()
        row3.addWidget(QLabel("Output prefix:"))
        self.sv_out = QLineEdit("vib")
        row3.addWidget(self.sv_out)
        vbox.addLayout(row3)

        row4 = QHBoxLayout()
        row4.addWidget(QLabel("Scale JSON (optional):"))
        scale_default = self._get_merlino_fit_root() / "scale_example.json"
        self.sv_scale = QLineEdit(str(scale_default) if scale_default.exists() else "")
        row4.addWidget(self.sv_scale)
        btn_scale = QPushButton("Browse")
        btn_scale.clicked.connect(self._browse_scale)
        row4.addWidget(btn_scale)
        vbox.addLayout(row4)

        row5 = QHBoxLayout()
        row5.addWidget(QLabel("Preview filter:"))
        self.sv_gic_filter = QComboBox()
        self.sv_gic_filter.addItems(["All", "Bonds", "Angles", "Dihedrals", "OOP", "Linear"])
        row5.addWidget(self.sv_gic_filter)
        vbox.addLayout(row5)

        self.run_survibfit_btn = QPushButton("Run Survibfit Vib")
        self.run_survibfit_btn.clicked.connect(lambda: AdvancedWindow.run_survibfit(self))
        vbox.addWidget(self.run_survibfit_btn)

        gic_group = QGroupBox("ReadGIC generator")
        gic_layout = QVBoxLayout(gic_group)
        rowg1 = QHBoxLayout()
        rowg1.addWidget(QLabel("XYZ file:"))
        self.sv_gic_xyz = QLineEdit(str(self.workdir / "xyzin"))
        rowg1.addWidget(self.sv_gic_xyz)
        btn_xyz = QPushButton("Browse")
        btn_xyz.clicked.connect(self._browse_gic_xyz)
        rowg1.addWidget(btn_xyz)
        gic_layout.addLayout(rowg1)

        rowg2 = QHBoxLayout()
        rowg2.addWidget(QLabel("Output file:"))
        self.sv_gic_out = QLineEdit(str(self.workdir / "gic.txt"))
        rowg2.addWidget(self.sv_gic_out)
        btn_out = QPushButton("Browse")
        btn_out.clicked.connect(self._browse_gic_out)
        rowg2.addWidget(btn_out)
        gic_layout.addLayout(rowg2)

        self.run_gic_btn = QPushButton("Generate ReadGIC")
        self.run_gic_btn.clicked.connect(lambda: AdvancedWindow.run_survibfit_gic(self))
        gic_layout.addWidget(self.run_gic_btn)

        self.preview_gic_btn = QPushButton("Preview ReadGIC")
        self.preview_gic_btn.clicked.connect(lambda: AdvancedWindow.preview_survibfit_gic(self))
        gic_layout.addWidget(self.preview_gic_btn)

        self.save_gic_btn = QPushButton("Save ReadGIC")
        self.save_gic_btn.clicked.connect(lambda: AdvancedWindow.save_survibfit_gic(self))
        gic_layout.addWidget(self.save_gic_btn)

        vbox.addWidget(gic_group)

        layout.addWidget(group)

    # ==================================================================
    # GF / PED and VPT2 / VCI panels
    # ==================================================================
    def _build_gf_shortcut_panel(self, layout: QVBoxLayout):
        group = QGroupBox("GF / PED")
        vbox = QVBoxLayout(group)

        info = QLabel(
            "Dedicated window for Wilson GF/PED from Cartesian Hessians, "
            "Merlino GIC definitions, B matrices, optional Pulay scaling, "
            "frequencies, normal modes and PED."
        )
        info.setWordWrap(True)
        vbox.addWidget(info)

        self.open_gf_btn = QPushButton("Open GF / PED window")
        self.open_gf_btn.clicked.connect(self.open_gf_window)
        vbox.addWidget(self.open_gf_btn)

        layout.addWidget(group)

    def open_gf_window(self):
        if not hasattr(self, "gf_window") or self.gf_window is None:
            self.gf_window = GFWindow(self.workdir, self.project_root, parent=self)
        self.gf_window.show()
        self.gf_window.raise_()
        self.gf_window.activateWindow()

    def _build_vpt2_vci_shortcut_panel(self, layout: QVBoxLayout):
        group = QGroupBox("VPT2 / VCI")
        vbox = QVBoxLayout(group)

        info = QLabel(
            "Dedicated window for anharmonic VPT2/VCI comparisons from "
            "canonical Merlino QFF inputs in Cartesian normal modes."
        )
        info.setWordWrap(True)
        vbox.addWidget(info)

        self.open_vpt2_vci_btn = QPushButton("Open VPT2 / VCI window")
        self.open_vpt2_vci_btn.clicked.connect(self.open_vpt2_vci_window)
        vbox.addWidget(self.open_vpt2_vci_btn)

        layout.addWidget(group)

    def open_vpt2_vci_window(self):
        if not hasattr(self, "vpt2_vci_window") or self.vpt2_vci_window is None:
            self.vpt2_vci_window = VPT2VCIWindow(self.workdir, self.project_root, parent=self)
        self.vpt2_vci_window.show()
        self.vpt2_vci_window.raise_()
        self.vpt2_vci_window.activateWindow()

    def _build_dvr_shortcut_panel(self, layout: QVBoxLayout):
        group = QGroupBox("Path DVR – Gaussian scan analysis")
        vbox = QVBoxLayout(group)

        info = QLabel(
            "DVR analysis has a dedicated window. It reads completed Gaussian "
            "scan/path logs and computes levels, profiles, properties, and "
            "optional Cremer-Pople labels."
        )
        info.setWordWrap(True)
        vbox.addWidget(info)

        self.open_dvr_btn = QPushButton("Open DVR window")
        self.open_dvr_btn.clicked.connect(self.open_dvr_window)
        vbox.addWidget(self.open_dvr_btn)

        layout.addWidget(group)

    def open_dvr_window(self):
        if not hasattr(self, "dvr_window") or self.dvr_window is None:
            self.dvr_window = DVRWindow(self.workdir, self.project_root, parent=self)
        self.dvr_window.show()
        self.dvr_window.raise_()
        self.dvr_window.activateWindow()

    # ==================================================================
    # Actions
    # ==================================================================

    def run_gicforge(self):
        charge, multiplicity = self._read_charge_multiplicity()

        writer = ProvinWriter(
            workdir=self.workdir,
            title=self.project_dir.name,
            charge=charge,
            multiplicity=multiplicity,
            keywords=self._collect_gicforge_keywords(),
        )
        writer.write()

        result = GICForgeLauncher(self.workdir).run()

        if not result.success:
            QMessageBox.critical(self, "GICForge failed", result.message)
            return

        gauin = self.workdir / "gauin"
        gauin_gjf = self.workdir / "gauin.gjf"
        if gauin.exists() and not gauin_gjf.exists():
            shutil.copyfile(gauin, gauin_gjf)
        self._ensure_gaussian_route_keywords_in_file(gauin)
        self._ensure_gaussian_route_keywords_in_file(gauin_gjf)

        self.run_msr_btn.setEnabled("msrin" in result.files)
        self.run_gaussian_btn.setEnabled("gauin" in result.files)

        self._export_project_files()

        self.prova_viewer = GICForgeViewer(
            self.workdir,
            run_gaussian_callback=self.run_gaussian,
            parent=self,
        )
        self.prova_viewer.show()

        QMessageBox.information(self, "Preparation OK", result.message)

    def run_prova(self):
        self.run_gicforge()

    def run_merlino_fit_prep(self):
        charge, multiplicity = self._read_charge_multiplicity()

        try:
            import sys
            merlino_fit_root = self._get_merlino_fit_root()
            if str(merlino_fit_root) not in sys.path:
                sys.path.insert(0, str(merlino_fit_root))
            from survibfit.modify_geom import read_xyz as read_xyz_local
            from survibfit.pipeline import primitives_from_topology, build_topology
            from survibfit.transforms import build_u_with_names, format_readgic_lines
            from survibfit.modify_geom import _load_topology_elements
        except Exception as e:
            QMessageBox.critical(self, "merlino_fit", f"Import error: {e}")
            return

        xyz_path = self.xyzin
        lines = self._generate_readgic_lines(
            read_xyz_local,
            _load_topology_elements,
            build_topology,
            primitives_from_topology,
            build_u_with_names,
            format_readgic_lines,
        )
        if lines is None:
            return
        lines = self._strip_readgic_values(lines)
        lines = self._ensure_full_stretches(
            lines,
            read_xyz_local,
            _load_topology_elements,
            primitives_from_topology,
        )
        lines = self._replace_ring_cyclic_coords(
            lines,
            read_xyz_local,
            _load_topology_elements,
            primitives_from_topology,
            build_topology,
        )

        route = self._build_gaussian_route_from_keywords()
        title = self.project_dir.name

        gauin = self.workdir / "gauin.gjf"
        with gauin.open("w") as fh:
            fh.write(f"%Nprocshared={self.config.gaussian_nproc}\n")
            fh.write(f"%Mem={self.config.gaussian_memory}\n")
            fh.write("%chk=gauin.chk\n")
            fh.write(route + "\n\n")
            fh.write(title + "\n\n")
            fh.write(f"{charge} {multiplicity}\n")

            atoms, coords_ang, _ = read_xyz_local(xyz_path)
            for sym, (x, y, z) in zip(atoms, coords_ang):
                fh.write(f"{sym:2s} {x: 11.6f} {y: 11.6f} {z: 11.6f}\n")
            fh.write("\n")
            for line in lines:
                fh.write(line + "\n")
            fh.write("\n")

        gauin_plain = self.workdir / "gauin"
        try:
            shutil.copyfile(gauin, gauin_plain)
            self._ensure_gaussian_route_keywords_in_file(gauin_plain)
        except Exception:
            pass
        self._ensure_gaussian_route_keywords_in_file(gauin)

        self.prova_viewer = GICForgeViewer(
            self.workdir,
            run_gaussian_callback=self.run_gaussian,
            parent=self,
        )
        self.prova_viewer.show()

        QMessageBox.information(self, "merlino_fit", f"Gaussian input written: {gauin}")

    def _strip_readgic_values(self, lines):
        import re
        out = []
        for line in lines:
            out.append(re.sub(r"\\(Value=[^)]+\\)", "", line))
        return out

    def _ensure_full_stretches(
        self,
        lines,
        read_xyz_local,
        load_topology_elements,
        primitives_from_topology,
    ):
        xyz_path = Path(self.xyzin)
        atoms, coords_ang, _ = read_xyz_local(xyz_path)
        coords_au = coords_ang / 0.52917721092
        atomic_number = load_topology_elements()
        Z = [atomic_number(a) for a in atoms]
        prims = primitives_from_topology(coords_au, Z, np.deg2rad(170.0))

        stretch_lines = []
        idx = 1
        for p in prims:
            if p.kind != "bond":
                continue
            i, j = p.atoms
            stretch_lines.append(f" Stre{idx:04d} =R({i+1:3d},{j+1:3d})")
            idx += 1

        rest = [ln for ln in lines if not ln.lstrip().startswith("Stre")]
        return stretch_lines + rest

    def _replace_ring_cyclic_coords(
        self,
        lines,
        read_xyz_local,
        load_topology_elements,
        primitives_from_topology,
        build_topology,
    ):
        try:
            from survibfit import transforms
            from survibfit.puckering_gaussian import ring_puckering_gic_lines
        except Exception:
            return lines

        xyz_path = Path(self.xyzin)
        atoms, coords_ang, _ = read_xyz_local(xyz_path)
        coords_au = coords_ang / 0.52917721092
        atomic_number = load_topology_elements()
        Z = [atomic_number(a) for a in atoms]

        prims = primitives_from_topology(coords_au, Z, np.deg2rad(170.0))
        _mol, _dg, ringset = build_topology(coords_au, Z)

        def _label_for_prim(p):
            idxs = [i + 1 for i in p.atoms]
            if p.kind == "angle":
                return f"A({idxs[0]:3d},{idxs[1]:3d},{idxs[2]:3d})"
            if p.kind == "dihedral":
                return f"D({idxs[0]:3d},{idxs[1]:3d},{idxs[2]:3d},{idxs[3]:3d})"
            return ""

        def _format_combo(prefix, U, idxs):
            out = []
            for col in range(U.shape[1]):
                terms = []
                for row, prim_idx in enumerate(idxs):
                    coeff = U[row, col]
                    if abs(coeff) < 1e-8:
                        continue
                    lab = _label_for_prim(prims[prim_idx])
                    if not lab:
                        continue
                    terms.append(f"{coeff: .5f}*{lab}")
                if not terms:
                    continue
                name = f"{prefix}{len(out)+1:04d}"
                expr = "[ " + " ".join(terms) + " ]"
                out.append(f" {name} ={expr}")
            return out

        U_ang, idx_ang = transforms._ring_cyclic_u(prims, ringset, "angle")
        U_dih, idx_dih = transforms._ring_cyclic_u(prims, ringset, "dihedral")
        U_cond, idx_cond = transforms.ring_condensed_dihedral_u(prims, coords_au, ringset, tol=1e-8, fd_step=1e-4)

        cvb_lines = _format_combo("CVB", U_ang, idx_ang)
        pucker_lines = []
        if ringset is not None:
            for ring in ringset:
                try:
                    if len(ring.atoms) >= 4:
                        pucker_lines.extend(ring_puckering_gic_lines(list(ring.atoms)))
                except Exception:
                    continue

        prefixes = ("CVB", "CTor", "T0", "RPck", "QPck", "PhiP")
        rest = [ln for ln in lines if not ln.lstrip().startswith(prefixes)]
        return rest + cvb_lines + pucker_lines

    def _ensure_gaussian_route_keywords_in_file(self, path: Path):
        """
        Ensure mandatory Gaussian route keywords exist in the first route line.
        """
        path = Path(path)
        if not path.exists():
            return
        try:
            lines = path.read_text(errors="replace").splitlines()
        except Exception:
            return

        route_idx = None
        for i, line in enumerate(lines):
            if line.strip().startswith("#"):
                route_idx = i
                break
        if route_idx is None:
            return

        route = lines[route_idx].strip()
        route_l = route.lower()
        changed = False
        if "iop(6/79=1,6/80=1)" not in route_l:
            route += " iop(6/79=1,6/80=1)"
            changed = True
        if "output=pickett" not in route_l:
            route += " output=Pickett"
            changed = True

        if changed:
            lines[route_idx] = route
            try:
                path.write_text("\n".join(lines) + "\n")
            except Exception:
                return

    def _build_gaussian_route_from_keywords(self) -> str:
        method_widget = self.kwd_widgets.get("METHOD")
        method = method_widget.currentText() if isinstance(method_widget, QComboBox) else ""

        loose = bool(self.kwd_widgets.get("LOOSE") and self.kwd_widgets["LOOSE"].isChecked())
        symmall = bool(self.kwd_widgets.get("SYMMALL") and self.kwd_widgets["SYMMALL"].isChecked())

        opt = bool(self.kwd_widgets.get("OPT") and self.kwd_widgets["OPT"].isChecked())
        harm = bool(self.kwd_widgets.get("HARM") and self.kwd_widgets["HARM"].isChecked())
        anh = bool(self.kwd_widgets.get("ANH") and self.kwd_widgets["ANH"].isChecked())

        if symmall:
            geom = "geom=(readallgic,gicallsymm)"
        else:
            geom = "geom=readallgic"

        parts = ["#P", geom]
        if loose:
            parts.append("Symm=loose")

        method_map = {
            "PCS0": "UFF",
            "PCS1": "HF3C",
            "HPCS2": "B3LYP EMPIRICALDISPERSION=GD4 6-31G*",
            "DPCS3": "revDSDPBEP86D4 gen",
        }
        if method in method_map:
            parts.append(method_map[method])

        if opt:
            parts.append("OPT")
        if harm:
            parts.append("FREQ")
        if anh:
            parts.append("ANHARM")

        route_l = " ".join(parts).lower()
        if "pop=cm5" not in route_l:
            parts.append("Pop=CM5")
        if "iop(6/79=1,6/80=1)" not in route_l:
            parts.append("iop(6/79=1,6/80=1)")
        if "output=pickett" not in route_l:
            parts.append("output=Pickett")

        return " ".join(parts)

    def run_msr(self):
        result = MSRLauncher(self.workdir).run()
        QMessageBox.information(self, "MSR", result.message)

    # ==============================================================
    # Gaussian (ASYNC, BACKGROUND)
    # ==============================================================

    def run_gaussian(self):
        self.gaussian_launcher = GaussianLauncher(self.workdir, parent=self, executable=self.config.gaussian_executable)
        self.gaussian_launcher.finished.connect(self._on_gaussian_finished)

        self.gaussian_launcher.start()

        QMessageBox.information(
            self,
            "Gaussian",
            "Gaussian started in background.\n"
            "You can continue using Merlino."
        )

    def _on_gaussian_finished(self, success: bool, message: str):
        if not success:
            QMessageBox.critical(self, "Gaussian failed", message)
            return

        self.gaussian_viewer = GaussianViewer(self.workdir)
        self.gaussian_viewer.show()

        QMessageBox.information(self, "Gaussian finished", message)
        self._export_project_files()



    def run_survibfit(self):
        log_path = Path(self.sv_log.text().strip())
        fchk_path = Path(self.sv_fchk.text().strip()) if self.sv_fchk.text().strip() else None
        scale_json = Path(self.sv_scale.text().strip()) if self.sv_scale.text().strip() else None
        out_prefix = self.sv_out.text().strip() or "vib"

        merlino_fit_root = self._get_merlino_fit_root()
        self.survibfit_launcher = SurvibfitLauncher(self.workdir, merlino_fit_root, parent=self)
        self.survibfit_launcher.finished.connect(self._on_survibfit_finished)
        self.survibfit_launcher.start_vib(log_path, fchk_path, out_prefix, scale_json)

    def run_survibfit_gic(self):
        xyz_path = Path(self.sv_gic_xyz.text().strip())
        out_path = Path(self.sv_gic_out.text().strip())
        merlino_fit_root = self._get_merlino_fit_root()
        self.survibfit_launcher = SurvibfitLauncher(self.workdir, merlino_fit_root, parent=self)
        self.survibfit_launcher.finished.connect(self._on_survibfit_finished)
        self.survibfit_launcher.start_gic(xyz_path, out_path, include_frag=False)

    def preview_survibfit_gic(self):
        from PySide6.QtWidgets import QDialog, QTextEdit, QVBoxLayout
        import sys
        merlino_fit_root = self._get_merlino_fit_root()
        if str(merlino_fit_root) not in sys.path:
            sys.path.insert(0, str(merlino_fit_root))
        from survibfit.modify_geom import read_xyz as read_xyz_local
        from survibfit.pipeline import _load_topology_elements, build_topology, primitives_from_topology
        from survibfit.transforms import build_u_with_names, format_readgic_lines

        xyz_path = Path(self.sv_gic_xyz.text().strip())
        if not xyz_path.exists():
            QMessageBox.critical(self, "Survibfit", f"XYZ not found: {xyz_path}")
            return

        lines = self._generate_readgic_lines(
            read_xyz_local, _load_topology_elements, build_topology,
            primitives_from_topology, build_u_with_names, format_readgic_lines
        )
        if lines is None:
            return
        filt = self.sv_gic_filter.currentText()
        if filt != "All":
            prefix = {
                "Bonds": "Stre",
                "Angles": "Bend",
                "Dihedrals": "Dihe",
                "OOP": "OuPl",
                "Linear": "LinB",
            }.get(filt, "")
            if prefix:
                lines = [l for l in lines if l.startswith(prefix)]

        dlg = QDialog(self)
        dlg.setWindowTitle("ReadGIC Preview")
        layout = QVBoxLayout(dlg)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText("\n".join(lines))
        layout.addWidget(text)
        dlg.resize(700, 500)
        dlg.exec()

    def save_survibfit_gic(self):
        import sys
        merlino_fit_root = self._get_merlino_fit_root()
        if str(merlino_fit_root) not in sys.path:
            sys.path.insert(0, str(merlino_fit_root))
        from survibfit.modify_geom import read_xyz as read_xyz_local
        from survibfit.pipeline import _load_topology_elements, build_topology, primitives_from_topology
        from survibfit.transforms import build_u_with_names, format_readgic_lines

        out_path = Path(self.sv_gic_out.text().strip())
        if not out_path.parent.exists():
            QMessageBox.critical(self, "Survibfit", f"Output folder not found: {out_path.parent}")
            return

        lines = self._generate_readgic_lines(
            read_xyz_local, _load_topology_elements, build_topology,
            primitives_from_topology, build_u_with_names, format_readgic_lines
        )
        if lines is None:
            return

        out_path.write_text("\n".join(lines) + "\n")
        QMessageBox.information(self, "Survibfit", f"ReadGIC saved to {out_path}")

    def _generate_readgic_lines(
        self,
        read_xyz_local,
        load_topology_elements,
        build_topology,
        primitives_from_topology,
        build_u_with_names,
        format_readgic_lines,
    ):
        xyz_path = Path(self.sv_gic_xyz.text().strip())
        if not xyz_path.exists():
            QMessageBox.critical(self, "Survibfit", f"XYZ not found: {xyz_path}")
            return None

        atoms, coords_ang, _ = read_xyz_local(xyz_path)
        coords_au = coords_ang / 0.52917721092
        atomic_number = load_topology_elements()
        Z = [atomic_number(a) for a in atoms]
        prims = primitives_from_topology(coords_au, Z, np.deg2rad(170.0))
        _, _, ringset = build_topology(coords_au, Z)
        _, names = build_u_with_names(prims, coords_au, Z=Z, ringset=ringset, include_frag=False)
        return format_readgic_lines(names)

    def _on_survibfit_finished(self, success: bool, message: str):
        if success:
            QMessageBox.information(self, "Survibfit", message)
        else:
            QMessageBox.critical(self, "Survibfit failed", message)

    def _browse_log(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Gaussian log", str(self.workdir), "Gaussian Logs (*.log *.out)"
        )
        if path:
            self.sv_log.setText(path)
            stem = Path(path).with_suffix("")
            for ext in (".fchk", ".fch"):
                cand = stem.with_suffix(ext)
                if cand.exists():
                    self.sv_fchk.setText(str(cand))
                    break

    def _browse_fchk(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Gaussian fchk", str(self.workdir), "Gaussian FCHK (*.fchk *.fch)"
        )
        if path:
            self.sv_fchk.setText(path)

    def _browse_scale(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select scale JSON", str(self.workdir), "JSON (*.json)"
        )
        if path:
            self.sv_scale.setText(path)

    def _browse_gic_xyz(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select XYZ", str(self.workdir), "XYZ (*.xyz *xyzin)"
        )
        if path:
            self.sv_gic_xyz.setText(path)

    def _browse_gic_out(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save ReadGIC", str(self.workdir / "gic.txt"), "Text (*.txt)"
        )
        if path:
            self.sv_gic_out.setText(path)

class AdvancedCalculationsWindow(AdvancedWindow):
    """Compatibility wrapper used by legacy tests/scripts."""

    def __init__(self, workdir: Path, project_dir: Path | None = None, parent=None):
        self._explicit_project_dir = Path(project_dir) if project_dir is not None else None
        super().__init__(workdir, parent=parent)

    def _init_project_dir(self):
        if self._explicit_project_dir is not None:
            self.project_name = self._explicit_project_dir.name
            self.project_out_dir = self._explicit_project_dir
            self.project_out_dir.mkdir(parents=True, exist_ok=True)
            return
        super()._init_project_dir()
