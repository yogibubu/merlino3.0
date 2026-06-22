from __future__ import annotations

import argparse
import configparser
import os
from pathlib import Path
import tempfile

import numpy as np

from .config import load_config, get_fit_config, get_u_config
from .terms import generate_terms, eval_terms, load_terms
from .fit import robust_fit
from .weights import scale_residuals
from .io import write_terms
from .pipeline import (
    primitives_from_topology,
    eval_primitives,
    b_matrix,
    load_u_matrix,
    q_from_s,
    gq_from_gx,
    g_matrix,
    g_matrix_derivs,
    build_topology,
)
from .transforms import build_u
from .vibrational_internal import modes_from_gaussian_log
from .modify_geom import read_xyz
from .pipeline import _load_topology_elements
from .transforms import build_u_with_names, format_readgic_lines
from .puckering_gaussian import (
    DEFAULT_GAUSSIAN_ROUTE,
    write_gaussian_scan_from_xyz,
)
from merlino_gic import define_gics_from_cartesian
from merlino_gic.model import parse_gicforge_line


def build_basis_cfg(cfg: configparser.ConfigParser, nvib: int):
    basis = {}
    section = cfg["basis"] if "basis" in cfg else {}
    default_mode = section.get("default_mode", "poly")
    for i in range(nvib):
        mode = section.get(f"mode_{i+1}", default_mode)
        params = {}
        if f"shift_{i+1}" in section:
            params["shift"] = float(section.get(f"shift_{i+1}"))
        if f"m_{i+1}" in section:
            params["m"] = int(section.get(f"m_{i+1}"))
        if f"a_{i+1}" in section:
            params["a"] = float(section.get(f"a_{i+1}"))
        if f"x0_{i+1}" in section:
            params["x0"] = float(section.get(f"x0_{i+1}"))
        basis[i] = {"mode": mode, "params": params}
    return basis


def term_indices(term, one_based=True):
    idxs = []
    for i, exp in term.exps:
        ii = i + 1 if one_based else i
        idxs.extend([ii] * exp)
    return idxs


def _fit_main(args):
    cfg = load_config(args.config)
    fit_cfg = get_fit_config(cfg)
    if args.cache_dir:
        import os
        os.environ["MERLINO_FIT_CACHE_DIR"] = args.cache_dir

    data = np.load(args.data)
    coords = data["coords"]
    Z = data["Z"]
    energy = data["energy"]
    grad = data["grad"]

    nat = coords.shape[1]
    linear_threshold = float(cfg.get("system", "linear_threshold", fallback=np.deg2rad(170.0)))

    prims = primitives_from_topology(coords[0], Z, linear_threshold)
    nprim = len(prims)

    u = None
    u_path = None
    u_mode = cfg.get("u", "mode", fallback="identity") if "u" in cfg else "identity"
    u_cfg = get_u_config(cfg)
    if "u" in cfg and "u_path" in cfg["u"]:
        u_path = cfg["u"]["u_path"]
    if u_path:
        u = load_u_matrix(u_path, nprim)
    elif u_mode == "auto":
        _, _, ringset = build_topology(coords[0], Z)
        def _group_filter_from_limit(spec):
            if spec is None:
                return None
            spec = str(spec).strip().lower()
            if spec in ("c", "cn", "cnv", "cnh"):
                return lambda lab: lab.startswith(("C", "sigma", "S", "E", "i"))
            if spec in ("d", "dn", "dnh", "dnd"):
                return lambda lab: lab.startswith(("C", "sigma", "S", "E", "i"))
            if spec == "poly":
                return lambda lab: lab.endswith(("_t", "_o", "_i")) or lab in ("E", "i")
            return None

        op_filter = _group_filter_from_limit(args.symmetry_group_limit or u_cfg.symmetry_group_limit)

        u = build_u(
            prims,
            coords[0],
            Z=Z,
            ringset=ringset,
            symmetry_mode=u_cfg.symmetry_mode,
            prune_mode=u_cfg.prune_mode,
            zeff_tol=u_cfg.zeff_tol,
            geometry_match_tol=u_cfg.geometry_match_tol,
            pattern_report_path=u_cfg.pattern_report_path,
            symmetrize_global=args.symmetrize_global or u_cfg.symmetrize_global,
            keep_a1_only=args.keep_a1_only or u_cfg.keep_a1_only,
            symmetry_tol=u_cfg.symmetry_tol,
            symmetry_max_n=u_cfg.symmetry_max_n,
            a1_tol=u_cfg.a1_tol,
            assign_symmetry_labels=args.assign_symmetry_labels or u_cfg.assign_symmetry_labels,
            symmetry_quasi_tol=args.symmetry_quasi_tol or u_cfg.symmetry_quasi_tol,
            symmetry_tol_h=args.symmetry_tol_h or u_cfg.symmetry_tol_H,
            heavy_only_orient=args.heavy_only_orient or u_cfg.heavy_only_orient,
            symmetry_center_idx=args.symmetry_center_idx if args.symmetry_center_idx is not None else u_cfg.symmetry_center_idx,
            ignore_isotopes=args.ignore_isotopes or u_cfg.ignore_isotopes,
            max_dev_strict=args.symmetry_max_dev_strict or u_cfg.symmetry_max_dev_strict,
            symmetry_tol_rel=args.symmetry_tol_rel or u_cfg.symmetry_tol_rel,
            symmetry_auto_max_n=args.symmetry_auto_max_n or u_cfg.symmetry_auto_max_n,
            symmetry_inertia_tol=args.symmetry_inertia_tol or u_cfg.symmetry_inertia_tol,
            symmetry_max_radius=args.symmetry_max_radius or u_cfg.symmetry_max_radius,
            symmetry_enforce_radial=not args.symmetry_no_radial_filter and u_cfg.symmetry_enforce_radial,
            symmetry_profile=args.symmetry_profile or u_cfg.symmetry_profile,
            symmetry_group_limit=op_filter,
            symmetry_confidence=args.symmetry_confidence or u_cfg.symmetry_confidence,
        )
    else:
        u = load_u_matrix(None, nprim)
    nvib = u.shape[1]

    basis_cfg = build_basis_cfg(cfg, nvib)
    term_list = None
    if "terms" in cfg and "term_list" in cfg["terms"]:
        term_list = cfg["terms"]["term_list"]
    terms = load_terms(term_list) if term_list else generate_terms(nvib)

    # Build design matrix
    rows = []
    targets = []

    for p in range(coords.shape[0]):
        s = eval_primitives(prims, coords[p])
        b = b_matrix(prims, coords[p], fit_cfg.fd_step)
        q = q_from_s(u, s)
        gx = grad[p].reshape(-1)
        gq = gq_from_gx(u, b, gx)

        phi, dphi = eval_terms(q, terms, basis_cfg)
        rows.append(phi)
        targets.append(energy[p])
        for i in range(nvib):
            rows.append(dphi[i, :])
            targets.append(gq[i])

    A = np.vstack(rows)
    y = np.array(targets)

    # Optional scaling pass (single-shot)
    scale_mode = cfg.get("fit", "scale", fallback="none")
    if scale_mode and scale_mode != "none":
        # initial residuals from unweighted fit
        coeff0 = robust_fit(
            A, y, delta=fit_cfg.delta, ridge=fit_cfg.ridge,
            max_iter=5, tol=fit_cfg.tol
        )
        res = y - A @ coeff0
        # split residuals into energy and gradient blocks
        nE = coords.shape[0]
        E_res = res[:nE]
        G_res = res[nE:]
        E_scaled, G_scaled, sE, sG = scale_residuals(E_res, G_res, mode=scale_mode)
        # scale A and y to balance energy/gradient blocks
        A[:nE, :] /= sE
        y[:nE] /= sE
        A[nE:, :] /= sG
        y[nE:] /= sG

    coeff = robust_fit(
        A,
        y,
        delta=fit_cfg.delta,
        ridge=fit_cfg.ridge,
        max_iter=fit_cfg.max_iter,
        tol=fit_cfg.tol,
    )

    # Write terms
    out_path = Path(args.out)
    write_terms(out_path, coeff, terms)


def _vib_main(args):
    scale_map = None
    if args.scale_json:
        import json
        scale_map = json.loads(Path(args.scale_json).read_text())
    elif args.scale_ini:
        import configparser
        cfg = configparser.ConfigParser()
        cfg.read(args.scale_ini)
        if "scale" in cfg:
            scale_map = {k: float(v) for k, v in cfg["scale"].items()}
        elif "u" in cfg:
            scale_map = {}
            for k, v in cfg["u"].items():
                if k.startswith("scale_"):
                    scale_map[k.replace("scale_", "")] = float(v)

    freqs, modes_q, U, prims = modes_from_gaussian_log(
        Path(args.log),
        fchk_path=Path(args.fchk) if args.fchk else None,
        scale_map=scale_map,
    )
    out = Path(args.out)
    np.save(out.with_suffix(".freqs.npy"), freqs)
    np.save(out.with_suffix(".modes_q.npy"), modes_q)
    np.save(out.with_suffix(".U.npy"), U)
    out.with_suffix(".freqs.txt").write_text("\n".join(f"{f: .6f}" for f in freqs) + "\n")


def _gic_main(args):
    if args.cache_dir:
        os.environ["MERLINO_FIT_CACHE_DIR"] = args.cache_dir
    if args.python_local and not _python_local_gic_allowed():
        raise SystemExit(
            "--python-local is a non-production diagnostic path. "
            "Set MERLINO_ALLOW_PYTHON_LOCAL_GIC=1 to use it explicitly."
        )
    atoms, coords_ang, _ = read_xyz(Path(args.xyz))
    if not args.python_local:
        symmetrize = bool(args.symmetrize_global or args.keep_a1_only or args.assign_symmetry_labels)
        workdir = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="survibfit_gicforge_"))
        definition = define_gics_from_cartesian(
            tuple(atoms),
            coords_ang,
            workdir=workdir,
            symmetrize=symmetrize,
        )
        lines = _canonical_gicforge_lines(definition.gaussian_input, definition.irreps, keep_a1_only=args.keep_a1_only)
        Path(args.out).write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    coords_au = coords_ang / 0.52917721092
    atomic_number = _load_topology_elements()
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)

    prims = primitives_from_topology(coords_au, Z, np.deg2rad(170.0))
    _, _, ringset = build_topology(coords_au, Z)

    def _group_filter_from_limit(spec):
        if spec is None:
            return None
        spec = str(spec).strip().lower()
        if spec in ("c", "cn", "cnv", "cnh"):
            return lambda lab: lab.startswith(("C", "sigma", "S", "E", "i"))
        if spec in ("d", "dn", "dnh", "dnd"):
            return lambda lab: lab.startswith(("C", "sigma", "S", "E", "i"))
        if spec == "poly":
            return lambda lab: lab.endswith(("_t", "_o", "_i")) or lab in ("E", "i")
        return None

    op_filter = _group_filter_from_limit(args.symmetry_group_limit)

    U, names = build_u_with_names(
        prims,
        coords_au,
        Z=Z,
        ringset=ringset,
        include_frag=args.include_frag,
        min_coeff=args.min_coeff,
        normalize=not args.no_normalize,
        symmetry_mode=args.symmetry_mode,
        prune_mode=args.prune_mode,
        zeff_tol=args.zeff_tol,
        geometry_match_tol=args.geometry_match_tol,
        pattern_report_path=args.pattern_report,
        symmetrize_global=args.symmetrize_global,
        keep_a1_only=args.keep_a1_only,
        assign_symmetry_labels=args.assign_symmetry_labels,
        symmetry_quasi_tol=args.symmetry_quasi_tol,
        symmetry_tol_h=args.symmetry_tol_h,
        heavy_only_orient=args.heavy_only_orient,
        symmetry_center_idx=args.symmetry_center_idx,
        ignore_isotopes=args.ignore_isotopes,
        max_dev_strict=args.symmetry_max_dev_strict,
        symmetry_tol_rel=args.symmetry_tol_rel,
        symmetry_auto_max_n=args.symmetry_auto_max_n,
        symmetry_inertia_tol=args.symmetry_inertia_tol,
        symmetry_max_radius=args.symmetry_max_radius,
        symmetry_enforce_radial=not args.symmetry_no_radial_filter,
        symmetry_profile=args.symmetry_profile,
        symmetry_group_limit=op_filter,
        symmetry_confidence=args.symmetry_confidence,
    )
    lines = format_readgic_lines(names)
    Path(args.out).write_text("\n".join(lines) + "\n")


def _canonical_gicforge_lines(text: str, irreps: tuple[str, ...], *, keep_a1_only: bool = False) -> list[str]:
    lines: list[str] = []
    irrep_index = 0
    for raw in text.splitlines():
        parsed = parse_gicforge_line(raw)
        if parsed is None:
            continue
        irrep = irreps[irrep_index] if irrep_index < len(irreps) else "UNK"
        irrep_index += 1
        if keep_a1_only and irrep not in {"A1", "A", "Ag", "A'"}:
            continue
        lines.append(raw.rstrip())
    return lines


def _python_local_gic_allowed() -> bool:
    value = os.environ.get("MERLINO_ALLOW_PYTHON_LOCAL_GIC", "")
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _pucker_gaussian_main(args):
    manifest = write_gaussian_scan_from_xyz(
        Path(args.xyz),
        args.ring,
        Path(args.gjf_out),
        phi_start=args.phi_start,
        phi_end=args.phi_end,
        phi_step=args.phi_step,
        charge=args.charge,
        multiplicity=args.multiplicity,
        mem=args.mem,
        nproc=args.nproc,
        chk_prefix=args.chk_prefix,
        title=args.title,
        route=args.route,
        constraint_mode=args.constraint_mode,
    )
    if args.manifest_out:
        import json
        Path(args.manifest_out).write_text(json.dumps(manifest, indent=2) + "\n")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_fit = sub.add_parser("fit")
    ap_fit.add_argument("--config", required=True)
    ap_fit.add_argument("--data", required=True)
    ap_fit.add_argument("--out", required=True)
    ap_fit.add_argument("--with-g", action="store_true")
    ap_fit.add_argument("--cache-dir", default=None, help="Enable disk cache for eval_primitives/b_matrix")
    ap_fit.add_argument("--symmetrize-global", action="store_true")
    ap_fit.add_argument("--keep-a1-only", action="store_true")
    ap_fit.add_argument("--assign-symmetry-labels", action="store_true")
    ap_fit.add_argument("--symmetry-quasi-tol", type=float, default=None)
    ap_fit.add_argument("--symmetry-tol-h", type=float, default=None)
    ap_fit.add_argument("--heavy-only-orient", action="store_true")
    ap_fit.add_argument("--symmetry-center-idx", type=int, default=None)
    ap_fit.add_argument("--ignore-isotopes", action="store_true")
    ap_fit.add_argument("--symmetry-max-dev-strict", type=float, default=None)
    ap_fit.add_argument("--symmetry-tol-rel", type=float, default=0.0)
    ap_fit.add_argument("--symmetry-auto-max-n", action="store_true")
    ap_fit.add_argument("--symmetry-inertia-tol", type=float, default=1.0e-3)
    ap_fit.add_argument("--symmetry-max-radius", type=float, default=None)
    ap_fit.add_argument("--symmetry-no-radial-filter", action="store_true")
    ap_fit.add_argument("--symmetry-profile", action="store_true")
    ap_fit.add_argument("--symmetry-group-limit", default=None)
    ap_fit.add_argument("--symmetry-confidence", action="store_true")

    ap_vib = sub.add_parser("vib")
    ap_vib.add_argument("--log", required=True)
    ap_vib.add_argument("--fchk", required=False)
    ap_vib.add_argument("--out", required=True)
    ap_vib.add_argument("--scale-json", required=False)
    ap_vib.add_argument("--scale-ini", required=False)

    ap_gic = sub.add_parser("gic")
    ap_gic.add_argument("--xyz", required=True)
    ap_gic.add_argument("--out", required=True)
    ap_gic.add_argument("--include-frag", action="store_true")
    ap_gic.add_argument("--min-coeff", type=float, default=1e-4)
    ap_gic.add_argument("--no-normalize", action="store_true")
    ap_gic.add_argument("--symmetry-mode", default="hybrid")
    ap_gic.add_argument("--prune-mode", default="svd")
    ap_gic.add_argument("--zeff-tol", type=float, default=0.05)
    ap_gic.add_argument("--geometry-match-tol", type=float, default=12.0)
    ap_gic.add_argument("--pattern-report", default=None)
    ap_gic.add_argument("--cache-dir", default=None, help="Enable disk cache for eval_primitives/b_matrix")
    ap_gic.add_argument("--workdir", default=None, help="GICForge work directory for canonical Python/Fortran output")
    ap_gic.add_argument("--python-local", action="store_true", help="Use the legacy pure-Python local GIC builder")
    ap_gic.add_argument("--symmetrize-global", action="store_true")
    ap_gic.add_argument("--keep-a1-only", action="store_true")
    ap_gic.add_argument("--assign-symmetry-labels", action="store_true")
    ap_gic.add_argument("--symmetry-quasi-tol", type=float, default=None)
    ap_gic.add_argument("--symmetry-tol-h", type=float, default=None)
    ap_gic.add_argument("--heavy-only-orient", action="store_true")
    ap_gic.add_argument("--symmetry-center-idx", type=int, default=None)
    ap_gic.add_argument("--ignore-isotopes", action="store_true")
    ap_gic.add_argument("--symmetry-max-dev-strict", type=float, default=None)
    ap_gic.add_argument("--symmetry-tol-rel", type=float, default=0.0)
    ap_gic.add_argument("--symmetry-auto-max-n", action="store_true")
    ap_gic.add_argument("--symmetry-inertia-tol", type=float, default=1.0e-3)
    ap_gic.add_argument("--symmetry-max-radius", type=float, default=None)
    ap_gic.add_argument("--symmetry-no-radial-filter", action="store_true")
    ap_gic.add_argument("--symmetry-profile", action="store_true")
    ap_gic.add_argument("--symmetry-group-limit", default=None)
    ap_gic.add_argument("--symmetry-confidence", action="store_true")

    ap_pucker = sub.add_parser("pucker-gaussian")
    ap_pucker.add_argument("--xyz", required=True)
    ap_pucker.add_argument("--ring", default=None)
    ap_pucker.add_argument("--gjf-out", required=True)
    ap_pucker.add_argument("--manifest-out", default=None)
    ap_pucker.add_argument("--phi-start", type=float, default=0.0)
    ap_pucker.add_argument("--phi-end", type=float, default=360.0)
    ap_pucker.add_argument("--phi-step", type=float, default=10.0)
    ap_pucker.add_argument(
        "--constraint-mode",
        choices=["functional-targets", "scan-to-zero"],
        default="functional-targets",
    )
    ap_pucker.add_argument("--charge", type=int, default=0)
    ap_pucker.add_argument("--multiplicity", type=int, default=1)
    ap_pucker.add_argument("--mem", default="16GB")
    ap_pucker.add_argument("--nproc", default="8")
    ap_pucker.add_argument("--chk-prefix", default="mw_path_phi")
    ap_pucker.add_argument("--title", default="Puckering constrained optimization")
    ap_pucker.add_argument("--route", default=DEFAULT_GAUSSIAN_ROUTE)

    args = ap.parse_args()
    if args.cmd == "fit":
        _fit_main(args)
    elif args.cmd == "vib":
        _vib_main(args)
    elif args.cmd == "gic":
        _gic_main(args)
    elif args.cmd == "pucker-gaussian":
        _pucker_gaussian_main(args)


if __name__ == "__main__":
    main()
