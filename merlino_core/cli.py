from __future__ import annotations

import argparse
import shutil
from pathlib import Path
import sys
import json

from merlino_core import build_run_manifest, ensure_workspace, load_config, write_default_config
from merlino_dvr import DVRRequest, build_path_analysis_args, write_dvr_manifest
from merlino_fortran.backends import BACKENDS, SOURCE_BACKENDS, resolve_backend, resolve_source_backend
from merlino_gaussian import summarize_gaussian_log
from merlino_gic import run_gicforge
from merlino_semiexp import (
    DEFAULT_SEMIEXP_OBSERVABLE,
    DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
    ParameterClassConstraint,
    QMParameterPredicate,
    SemiexperimentalFitRequest,
    fit_semiexperimental_geometry,
    read_observations,
    semiexperimental_latex_tables,
    write_semiexperimental_html_report,
)
from merlino_vpt2_vci import (
    QuarticForceField,
    VCIOptions,
    load_force_field,
    run_gf_report_from_fchk,
    run_vpt2_vci_report,
    solve_vci,
    write_csv_tables,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="merlino", description="Merlino4 workflow CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    init = sub.add_parser("init", help="Create canonical workspace folders and merlino.toml")
    init.add_argument("workdir", type=Path)
    init.add_argument("--overwrite-config", action="store_true")

    config = sub.add_parser("config", help="Show resolved configuration")
    config.add_argument("--workdir", type=Path, default=Path("."))
    config.add_argument("--config", type=Path)

    gf = sub.add_parser("gf", help="Run GF/PED from an FCHK adapter")
    gf.add_argument("--fchk", type=Path, required=True)
    gf.add_argument("--out", type=Path)
    gf.add_argument("--run-dir", type=Path)
    gf.add_argument("--csv-dir", type=Path)

    vci = sub.add_parser("vci", help="Run VPT2/VCI from canonical QFF and optional FCHK frequencies")
    vci.add_argument("--qff", type=Path)
    vci.add_argument("--fchk", type=Path)
    vci.add_argument("--max-quanta", type=int, default=2)
    vci.add_argument("--roots", type=int, default=6)
    vci.add_argument("--active-modes", default="")
    vci.add_argument("--out", type=Path)
    vci.add_argument("--run-dir", type=Path)
    vci.add_argument("--csv-dir", type=Path)

    gic = sub.add_parser("gic", help="Run GICForge in a work directory")
    gic.add_argument("--workdir", type=Path, required=True)
    gic.add_argument("--executable", type=Path)

    summary = sub.add_parser("gaussian-summary", help="Summarize a Gaussian log/out file")
    summary.add_argument("log", type=Path)

    sub.add_parser("backends", help="Show configured backend availability")
    sub.add_parser("compare-backends", help="Run small backend consistency checks")

    dvr = sub.add_parser("dvr-args", help="Build DVR command args and manifest without executing")
    dvr.add_argument("--repo-root", type=Path, required=True)
    dvr.add_argument("--log", type=Path, required=True)
    dvr.add_argument("--outdir", type=Path, required=True)
    dvr.add_argument("--figdir", type=Path, required=True)
    dvr.add_argument("--prefix", default="puckering_dvr")
    dvr.add_argument("--boundary", default="periodic")
    dvr.add_argument("--solver", default="fourier")
    dvr.add_argument("--no-rotconst", action="store_true")
    dvr.add_argument("--label-cremer-pople", action="store_true")

    semiexp = sub.add_parser(
        "semiexp",
        help="Fit semiexperimental equilibrium geometry with the Cartesian/GIC Merlino standard solver",
    )
    semiexp.add_argument("--xyz", type=Path, required=True, help="Initial parent Cartesian geometry in XYZ format")
    semiexp.add_argument("--observations", type=Path, required=True, help="CSV/JSON/TOML with isotopologue B0 constants and corrections")
    semiexp.add_argument("--outdir", type=Path, required=True, help="Output directory for geometry, parameters, residuals and manifest")
    semiexp.add_argument("--backend", choices=("python", "fortran77"), default="python", help="Numerical backend requested by CLI/GUI")
    semiexp.add_argument("--fixed", default="", help="Comma/semicolon-separated GIC label substrings to keep fixed")
    semiexp.add_argument("--max-iter", type=int, default=20, help="Maximum LM iterations; default is capped for semiexp fits")
    semiexp.add_argument("--step", type=float, default=1.0e-4, help="Finite step for rotational observable derivatives with respect to GICs")
    semiexp.add_argument("--damping", type=float, default=1.0e-8, help="Initial Levenberg-Marquardt damping")
    semiexp.add_argument("--max-step", type=float, default=0.25, help="Maximum active-GIC step norm per iteration")
    semiexp.add_argument(
        "--observable",
        choices=("moments", "rotational_constants", "auto"),
        default=DEFAULT_SEMIEXP_OBSERVABLE,
        help="Fit target; default moments is the Merlino standard because it is more stable than reciprocal rotational constants",
    )
    semiexp.add_argument(
        "--rotational-components",
        choices=("auto", "ABC", "AB", "AC", "BC"),
        default=DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
        help="Rotational constants to use when observable=rotational_constants; auto chooses best-conditioned pair for planar molecules",
    )
    semiexp.add_argument(
        "--qm-predicate",
        action="append",
        default=[],
        help="QM prior as label_pattern:value:sigma[:source]; can be repeated",
    )
    semiexp.add_argument(
        "--parameter-class",
        action="append",
        default=[],
        help="Class constraint as name:shared|fixed:pattern[|pattern...]; can be repeated",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "init":
        layout = ensure_workspace(args.workdir)
        config_path = write_default_config(layout.root / "merlino.toml", overwrite=args.overwrite_config)
        print(f"workspace: {layout.root}")
        print(f"config: {config_path}")
        return 0

    if args.command == "config":
        config = load_config(args.config, workdir=args.workdir)
        for key, value in config.to_dict().items():
            print(f"{key}: {value}")
        return 0

    if args.command == "gf":
        report = run_gf_report_from_fchk(args.fchk)
        out = args.out or Path("gf_ped_report.txt")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report.text + "\n", encoding="utf-8")
        outputs = {"report": out}
        if args.csv_dir is not None:
            outputs.update({f"csv_{name}": path for name, path in write_csv_tables(report, args.csv_dir).items()})
        run_dir = args.run_dir or out.parent
        build_run_manifest(
            workflow="gf",
            status="completed",
            run_dir=run_dir,
            inputs={"fchk": args.fchk},
            outputs=outputs,
            backend={"adapter": "gaussian-fchk", "solver": "python"},
        ).write(Path(run_dir) / "gf_manifest.json")
        print(out)
        return 0

    if args.command == "vci":
        active_modes = _parse_active_modes(args.active_modes)
        qff = load_force_field(fchk_path=args.fchk, qff_path=args.qff)
        report = run_vpt2_vci_report(
            qff,
            max_quanta=args.max_quanta,
            roots=args.roots,
            options=VCIOptions(active_modes=active_modes),
        )
        out = args.out or Path("vpt2_vci_report.txt")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report.text + "\n", encoding="utf-8")
        outputs = {"report": out}
        if args.csv_dir is not None:
            outputs.update({f"csv_{name}": path for name, path in write_csv_tables(report, args.csv_dir).items()})
        run_dir = args.run_dir or out.parent
        inputs = {}
        if args.fchk is not None:
            inputs["fchk"] = args.fchk
        if args.qff is not None:
            inputs["qff"] = args.qff
        build_run_manifest(
            workflow="vpt2_vci",
            status="completed",
            run_dir=run_dir,
            inputs=inputs,
            outputs=outputs,
            parameters={"max_quanta": args.max_quanta, "roots": args.roots, "active_modes": active_modes},
            backend={"solver": "python"},
        ).write(Path(run_dir) / "vpt2_vci_manifest.json")
        print(out)
        return 0

    if args.command == "dvr-args":
        request = DVRRequest(
            repo_root=args.repo_root,
            log_path=args.log,
            outdir=args.outdir,
            figdir=args.figdir,
            prefix=args.prefix,
            boundary=args.boundary,
            solver=args.solver,
            compute_rotconst=not args.no_rotconst,
            label_cremer_pople=args.label_cremer_pople,
        )
        request.outdir.mkdir(parents=True, exist_ok=True)
        request.figdir.mkdir(parents=True, exist_ok=True)
        dvr_args = build_path_analysis_args(request)
        manifest = write_dvr_manifest(request, dvr_args)
        print(" ".join(dvr_args))
        print(f"manifest: {manifest}")
        return 0

    if args.command == "gic":
        result = run_gicforge(args.workdir, executable=args.executable)
        print(f"manifest: {result.manifest}")
        for name, path in sorted(result.files.items()):
            print(f"{name}: {path}")
        return 0

    if args.command == "semiexp":
        fixed = _parse_fixed_parameters(args.fixed)
        observations = read_observations(args.observations)
        request = SemiexperimentalFitRequest(
            initial_geometry=args.xyz,
            observations=observations,
            fixed_parameters=fixed,
            observable=args.observable,
            rotational_components=args.rotational_components,
            qm_predicates=_parse_qm_predicates(args.qm_predicate),
            parameter_classes=_parse_parameter_classes(args.parameter_class),
        )
        result = fit_semiexperimental_geometry(
            request,
            max_iter=args.max_iter,
            step=args.step,
            damping=args.damping,
            max_step=args.max_step,
            outdir=args.outdir,
        )
        report_path = write_semiexperimental_html_report(args.outdir / "semiexp_report.html", result, request)
        tables_path = args.outdir / "semiexp_tables.tex"
        tables = semiexperimental_latex_tables(result)
        tables_path.write_text(
            "\n\n".join(f"% {name}\n{table}" for name, table in tables.items()),
            encoding="utf-8",
        )
        _append_manifest_output(args.outdir / "semiexp_manifest.json", "html_report", report_path)
        _append_manifest_output(args.outdir / "semiexp_manifest.json", "latex_tables", tables_path)
        print(f"manifest: {result.manifest}")
        print(f"report: {report_path}")
        print(f"rms_MHz: {result.rms_MHz:.8g}")
        print(f"iterations: {result.iterations}")
        print(f"stationary_point: {result.stationary_point}")
        print(f"convergence: {result.diagnostics.convergence_reason}")
        print(f"rank: {result.diagnostics.rank}")
        print(f"condition_number: {result.diagnostics.condition_number:.8g}")
        print(f"observable: {result.diagnostics.observable}")
        print(f"components: {','.join(result.diagnostics.components)}")
        print(f"backend: {args.backend}")
        return 0

    if args.command == "gaussian-summary":
        summary = summarize_gaussian_log(args.log)
        print(f"path: {summary.path}")
        print(f"normal_termination: {summary.normal_termination}")
        print(f"scf_count: {len(summary.scf_energies_hartree)}")
        if summary.scf_energies_hartree:
            print(f"last_scf_hartree: {summary.scf_energies_hartree[-1]}")
        print(f"standard_orientation_count: {summary.standard_orientation_count}")
        print(f"input_orientation_count: {summary.input_orientation_count}")
        print(f"scan_marker_count: {summary.scan_marker_count}")
        print(f"puckering_marker_count: {summary.puckering_marker_count}")
        print(f"frequency_count: {len(summary.frequencies_cm)}")
        print(f"last_orientation_atoms: {len(summary.last_orientation)}")
        return 0

    if args.command == "backends":
        for name in sorted(BACKENDS):
            try:
                path = resolve_backend(name)
                status = "available"
            except Exception as exc:
                path = exc
                status = "missing"
            print(f"{name}: {status} ({path})")
        for name in sorted(SOURCE_BACKENDS):
            try:
                path = resolve_source_backend(name)
                status = "available"
            except Exception as exc:
                path = exc
                status = "missing"
            print(f"{name}: {status} ({path})")
        gaussian = shutil.which(load_config().gaussian_executable)
        print(f"gaussian: {'available' if gaussian else 'missing'} ({gaussian or load_config().gaussian_executable})")
        return 0

    if args.command == "compare-backends":
        import numpy as np

        qff = QuarticForceField(
            harmonic_frequencies_cm=np.array([1000.0, 1500.0]),
            cubic_cm={(0, 0, 1): -2.0},
            quartic_cm={(0, 0, 0, 0): 0.8},
        )
        dense = solve_vci(qff, max_quanta=3, n_roots=3, method="dense")
        davidson = solve_vci(qff, max_quanta=3, n_roots=3, method="davidson")
        delta = float(np.max(np.abs(dense.energies_cm - davidson.energies_cm)))
        print(f"python_dense_vs_davidson_max_delta_cm-1: {delta:.6g}")
        try:
            source = resolve_source_backend("vpt2_vci")
            print(f"fortran_vpt2_vci_source: available ({source})")
        except Exception as exc:
            print(f"fortran_vpt2_vci_source: missing ({exc})")
        return 0 if delta < 1.0e-6 else 1

    return 2


def _parse_active_modes(raw: str) -> tuple[int, ...] | None:
    if not raw.strip():
        return None
    values = tuple(int(part.strip()) - 1 for part in raw.replace(";", ",").split(",") if part.strip())
    if any(value < 0 for value in values):
        raise ValueError("active modes are one-based")
    return values


def _parse_fixed_parameters(raw: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in raw.replace(";", ",").split(",") if part.strip())


def _parse_qm_predicates(items: list[str]) -> tuple[QMParameterPredicate, ...]:
    predicates = []
    for item in items:
        parts = item.split(":")
        if len(parts) not in {3, 4}:
            raise ValueError("--qm-predicate must be label_pattern:value:sigma[:source]")
        source = parts[3] if len(parts) == 4 else "qm"
        predicates.append(QMParameterPredicate(parts[0], float(parts[1]), float(parts[2]), source=source))
    return tuple(predicates)


def _parse_parameter_classes(items: list[str]) -> tuple[ParameterClassConstraint, ...]:
    constraints = []
    for item in items:
        parts = item.split(":", 2)
        if len(parts) != 3:
            raise ValueError("--parameter-class must be name:shared|fixed:pattern[|pattern...]")
        patterns = tuple(part.strip() for part in parts[2].split("|") if part.strip())
        constraints.append(ParameterClassConstraint(parts[0].strip(), patterns, parts[1].strip()))
    return tuple(constraints)


def _append_manifest_output(manifest_path: Path, name: str, path: Path) -> None:
    if not manifest_path.exists():
        return
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    data.setdefault("outputs", {})[name] = str(path)
    if path.is_file():
        from merlino_core import sha256_file

        data.setdefault("output_sha256", {})[name] = sha256_file(path)
    manifest_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
