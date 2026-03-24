from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QComboBox, QDialog, QDialogButtonBox, QFormLayout, QMessageBox

from geometry.thermo_trasl import read_xyz_from_xyzin
from merlino_fit.survibfit.modify_geom import (
    ANG_TO_BOHR,
    BOHR_TO_ANG,
    _load_topology_elements,
    _disable_bdpcs3_fit,
    bdpcs3_delta_and_order,
    bdpcs3_delta_and_order_updated,
    topology_bond_order_for_pair,
    primitives_from_topology,
    eval_primitives,
    isotopic_masses_au,
    _backtransform_iterative,
    rotational_constants,
    write_xyz,
)


@dataclass
class Bdpcs3Result:
    symbols: list[str]
    coords_dpcs3_ang: np.ndarray
    coords_bdpcs3_ang: np.ndarray
    bond_orders: list[float]
    bond_targets: list[tuple[int, int, float, float]]  # i, j, r_corr(Ang), bond order
    s0_bond: np.ndarray
    s1_bond: np.ndarray
    bond_target_lengths: list[float]
    bond_residuals: list[float]
    backtransform_resid_norm: float
    backtransform_max_abs: float
    worst_primitives: list[tuple[str, tuple[int, ...], float]]
    backtransform_used_weights: bool
    B0_mhz: np.ndarray
    B1_mhz: np.ndarray


def ask_bdpcs3_version(parent, current_version: str) -> str | None:
    dlg = QDialog(parent)
    dlg.setWindowTitle("BDPCS3 settings")
    layout = QFormLayout(dlg)
    combo = QComboBox(dlg)
    combo.addItems(["legacy", "updated"])
    combo.setCurrentText(current_version)
    layout.addRow("BDPCS3 version", combo)

    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dlg.accept)
    buttons.rejected.connect(dlg.reject)
    layout.addRow(buttons)

    if dlg.exec() != QDialog.Accepted:
        return None
    return combo.currentText()


def write_xyzin_with_coords(src_xyzin_path: Path, out_path: Path, symbols, coords_ang):
    lines = []
    if src_xyzin_path.exists():
        src_lines = src_xyzin_path.read_text().splitlines()
    else:
        src_lines = []

    nat = len(symbols)
    comment = "BDPCS3"
    if len(src_lines) >= 2:
        try:
            int(src_lines[0].strip())
            comment = src_lines[1]
        except Exception:
            pass

    lines.append(str(nat))
    lines.append(comment)
    for sym, (x, y, z) in zip(symbols, coords_ang):
        lines.append(f"{sym} {x: .8f} {y: .8f} {z: .8f}")

    i = 0
    if src_lines:
        try:
            n0 = int(src_lines[0])
            i = n0 + 2
        except Exception:
            i = 0
    lines.extend(src_lines[i:])
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def compute_bdpcs3(
    xyzin_path: Path,
    bdpcs3_version: str,
) -> Bdpcs3Result:
    symbols, coords = read_xyz_from_xyzin(str(xyzin_path))
    if not symbols:
        raise ValueError("No atoms found in xyzin.")

    coords_ang = np.array(coords, dtype=float)
    coords_au = coords_ang * ANG_TO_BOHR

    atomic_number = _load_topology_elements()
    Z = np.array([atomic_number(sym) for sym in symbols], dtype=int)

    masses, _ = isotopic_masses_au(Z)

    prims_all = primitives_from_topology(coords_au, Z, linear_threshold=np.deg2rad(170.0))
    prims_bond = [p for p in prims_all if p.kind == "bond"]
    if not prims_bond:
        raise ValueError("No bond primitives detected.")

    s0_bond = eval_primitives(prims_bond, coords_au)
    s_all = eval_primitives(prims_all, coords_au)
    s_target_all = s_all.copy()

    _disable_bdpcs3_fit()
    bdpcs3_fn = bdpcs3_delta_and_order_updated if bdpcs3_version == "updated" else bdpcs3_delta_and_order

    bond_orders = []
    bond_targets = []
    bond_target_lengths = []
    bo_cache = {}
    for idx, p in enumerate(prims_all):
        if p.kind != "bond":
            continue
        i, j = p.atoms
        r_ang = s_all[idx] * BOHR_TO_ANG
        bo_topo = topology_bond_order_for_pair(i, j, Z, coords_ang, cache=bo_cache)
        dlt_r, bndord = bdpcs3_fn(
            int(Z[i]),
            int(Z[j]),
            r_ang,
            bond_order_override=bo_topo,
        )
        r_corr = r_ang + dlt_r
        s_target_all[idx] = r_corr * ANG_TO_BOHR
        bond_targets.append((i, j, r_corr, bndord))
        bond_orders.append(bndord)
        bond_target_lengths.append(r_corr)

    def _eval_backtransform(coords_bdpcs3_au):
        s1_bond = eval_primitives(prims_bond, coords_bdpcs3_au)
        coords_bdpcs3_ang = coords_bdpcs3_au * BOHR_TO_ANG
        s1_all = eval_primitives(prims_all, coords_bdpcs3_au)
        resid = s_target_all - s1_all
        backtransform_resid_norm = float(np.linalg.norm(resid))
        backtransform_max_abs = float(np.max(np.abs(resid)))
        worst_primitives = []
        if resid.size:
            idxs = np.argsort(np.abs(resid))[-8:][::-1]
            for idx in idxs:
                p = prims_all[idx]
                diff = float(resid[idx])
                if p.kind == "bond":
                    diff = diff * BOHR_TO_ANG
                worst_primitives.append((p.kind, tuple(p.atoms), diff))
        bond_residuals = []
        max_bond_resid = 0.0
        for r1, r_target in zip(s1_bond, bond_target_lengths):
            resid_bond = (r1 * BOHR_TO_ANG) - r_target
            bond_residuals.append(resid_bond)
            if abs(resid_bond) > max_bond_resid:
                max_bond_resid = abs(resid_bond)
        return (
            s1_bond,
            coords_bdpcs3_ang,
            bond_residuals,
            backtransform_resid_norm,
            backtransform_max_abs,
            worst_primitives,
            max_bond_resid,
        )

    coords_bdpcs3_au = _backtransform_iterative(
        s_target_all,
        coords_au,
        prims_all,
        masses=masses,
        max_iter=50,
        tol=1e-8,
        damping=1.0,
        adaptive=True,
    )
    (
        s1_bond,
        coords_bdpcs3_ang,
        bond_residuals,
        backtransform_resid_norm,
        backtransform_max_abs,
        worst_primitives,
        max_bond_resid,
    ) = _eval_backtransform(coords_bdpcs3_au)

    used_weights = False
    if max_bond_resid > 0.02:
        weights = [1.0 if p.kind == "bond" else 0.2 for p in prims_all]
        coords_bdpcs3_au = _backtransform_iterative(
            s_target_all,
            coords_au,
            prims_all,
            masses=masses,
            max_iter=80,
            tol=1e-8,
            damping=0.8,
            adaptive=True,
            weights=weights,
        )
        (
            s1_bond,
            coords_bdpcs3_ang,
            bond_residuals,
            backtransform_resid_norm,
            backtransform_max_abs,
            worst_primitives,
            _,
        ) = _eval_backtransform(coords_bdpcs3_au)
        used_weights = True

    B0 = rotational_constants(coords_au, masses) * 1000.0
    B1 = rotational_constants(coords_bdpcs3_au, masses) * 1000.0

    return Bdpcs3Result(
        symbols=symbols,
        coords_dpcs3_ang=coords_ang,
        coords_bdpcs3_ang=coords_bdpcs3_ang,
        bond_orders=bond_orders,
        bond_targets=bond_targets,
        s0_bond=s0_bond,
        s1_bond=s1_bond,
        bond_target_lengths=bond_target_lengths,
        bond_residuals=bond_residuals,
        backtransform_resid_norm=backtransform_resid_norm,
        backtransform_max_abs=backtransform_max_abs,
        worst_primitives=worst_primitives,
        backtransform_used_weights=used_weights,
        B0_mhz=B0,
        B1_mhz=B1,
    )


def write_bdpcs3_outputs(
    working_dir: Path,
    xyzin_path: Path,
    result: Bdpcs3Result,
    bdpcs3_version: str,
    basic_charge: int,
    basic_mult: int,
):
    report_path = working_dir / "bdpcs3.report"
    bdpcs3_xyz_path = working_dir / "bdpcs3.xyz"
    bdpcs3_xyzin_path = working_dir / "bdpcs3.xyzin"
    bdpcs3_gjf_path = working_dir / "BDPCS3.gjf"

    lines = []
    lines.append("DPCS3 to BDPCS3 report")
    lines.append(f"Source: {xyzin_path}")
    lines.append(f"BDPCS3 version: {bdpcs3_version}")
    lines.append("Units: Angstrom, MHz")
    lines.append("")
    lines.append("Bond lengths (Angstrom)")
    lines.append("")
    lines.append("  i  El |  j  El |        DPCS3 |      Target |       BDPCS3 |     Delta |  Resid | BondOrder")
    lines.append("--------+--------+-------------+-------------+-------------+----------+--------+----------")
    for idx, (r0, r1) in enumerate(zip(result.s0_bond, result.s1_bond)):
        i, j = result.bond_targets[idx][0], result.bond_targets[idx][1]
        delta = (r1 - r0) * BOHR_TO_ANG
        target = result.bond_target_lengths[idx]
        resid = result.bond_residuals[idx]
        lines.append(
            f"{i+1:3d} {result.symbols[i]:2s} | {j+1:3d} {result.symbols[j]:2s} |"
            f" {r0 * BOHR_TO_ANG: 11.6f} | {target: 11.6f} | {r1 * BOHR_TO_ANG: 11.6f} |"
            f" {delta: 8.6f} | {resid: 6.3f} | {result.bond_orders[idx]:8.3f}"
        )
    lines.append("")
    lines.append("Backtransform residuals (internal units; bonds in Angstrom)")
    lines.append(f"  Weighted backtransform: {'yes' if result.backtransform_used_weights else 'no'}")
    lines.append(f"  L2 norm: {result.backtransform_resid_norm:.6f}")
    lines.append(f"  Max abs: {result.backtransform_max_abs:.6f}")
    if result.worst_primitives:
        lines.append("  Worst primitives (kind, atoms, residual):")
        for kind, atoms, diff in result.worst_primitives:
            if kind == "bond":
                lines.append(f"    {kind} {tuple(a + 1 for a in atoms)} {diff: .6f} A")
            else:
                lines.append(f"    {kind} {tuple(a + 1 for a in atoms)} {diff: .6f}")
    lines.append("")
    lines.append("Rotational constants (MHz)")
    lines.append("")
    lines.append("         |           A |           B |           C")
    lines.append("---------+-------------+-------------+-------------")
    lines.append(f"DPCS3    | {result.B0_mhz[0]: 11.6f} | {result.B0_mhz[1]: 11.6f} | {result.B0_mhz[2]: 11.6f}")
    lines.append(f"BDPCS3   | {result.B1_mhz[0]: 11.6f} | {result.B1_mhz[1]: 11.6f} | {result.B1_mhz[2]: 11.6f}")
    lines.append("---------+-------------+-------------+-------------")
    lines.append("")
    lines.append("DPCS3 Cartesian coordinates (Angstrom)")
    lines.append(f"{len(result.symbols)}")
    lines.append("  #  El |           X |           Y |           Z")
    lines.append("------+----+-------------+-------------+-------------")
    for idx, (sym, (x, y, z)) in enumerate(zip(result.symbols, result.coords_dpcs3_ang), start=1):
        lines.append(f"{idx:3d} {sym:2s} | {x: 11.6f} | {y: 11.6f} | {z: 11.6f}")
    lines.append("")
    lines.append("BDPCS3 Cartesian coordinates (Angstrom)")
    lines.append(f"{len(result.symbols)}")
    lines.append("  #  El |           X |           Y |           Z")
    lines.append("------+----+-------------+-------------+-------------")
    for idx, (sym, (x, y, z)) in enumerate(zip(result.symbols, result.coords_bdpcs3_ang), start=1):
        lines.append(f"{idx:3d} {sym:2s} | {x: 11.6f} | {y: 11.6f} | {z: 11.6f}")
    lines.append("")
    lines.append(f"BDPCS3 XYZ: {bdpcs3_xyz_path}")
    lines.append(f"BDPCS3 xyzin: {bdpcs3_xyzin_path}")
    lines.append(f"Gaussian input: {bdpcs3_gjf_path}")
    lines.append("")

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_xyz(bdpcs3_xyz_path, result.symbols, result.coords_bdpcs3_ang, comment="BDPCS3")
    write_xyzin_with_coords(xyzin_path, bdpcs3_xyzin_path, result.symbols, result.coords_bdpcs3_ang)

    gjf_lines = []
    gjf_lines.append("#HF geom=modredundant Pop=CM5 iop(6/79=1,6/80=1) output=Pickett")
    gjf_lines.append("")
    gjf_lines.append("DPCS3 coords + BDPCS3 bond constraints")
    gjf_lines.append("")
    gjf_lines.append(f"{basic_charge} {basic_mult}")
    for sym, (x, y, z) in zip(result.symbols, result.coords_dpcs3_ang):
        gjf_lines.append(f"{sym:2s} {x: 11.6f} {y: 11.6f} {z: 11.6f}")
    gjf_lines.append("")
    for (i, j, _, _), r_target in zip(result.bond_targets, result.bond_target_lengths):
        gjf_lines.append(f"B {i+1} {j+1} {r_target:.6f}")
    gjf_lines.append("")
    bdpcs3_gjf_path.write_text("\n".join(gjf_lines) + "\n", encoding="utf-8")

    return report_path, bdpcs3_xyz_path, bdpcs3_xyzin_path, bdpcs3_gjf_path


def run_bdpcs3_report(
    parent,
    working_dir: Path,
    xyzin_path: Path,
    basic_charge: int,
    basic_mult: int,
    bdpcs3_version: str,
):
    if not xyzin_path.exists():
        QMessageBox.warning(parent, "BDPCS3", "xyzin not found in working.")
        return None

    version = ask_bdpcs3_version(parent, bdpcs3_version)
    if version is None:
        return None

    try:
        result = compute_bdpcs3(xyzin_path, version)
    except Exception as e:
        QMessageBox.critical(parent, "BDPCS3", str(e))
        return None

    report_path = xyz_path = xyzin_out_path = gjf_path = None
    try:
        report_path, xyz_path, xyzin_out_path, gjf_path = write_bdpcs3_outputs(
            working_dir,
            xyzin_path,
            result,
            version,
            basic_charge,
            basic_mult,
        )
    except Exception as e:
        QMessageBox.warning(parent, "BDPCS3", f"Report written with warnings: {e}")

    return {
        "version": version,
        "report_path": report_path,
        "xyz_path": xyz_path,
        "xyzin_path": xyzin_out_path,
        "gjf_path": gjf_path,
    }
