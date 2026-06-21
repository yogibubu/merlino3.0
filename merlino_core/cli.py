from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path
import sys
import json

from merlino_core import build_run_manifest, ensure_workspace, load_config, write_default_config
from merlino_dvr import DVRRequest, build_path_analysis_args, write_dvr_manifest
from merlino_fortran.backends import BACKENDS, SOURCE_BACKENDS, resolve_backend, resolve_source_backend
from merlino_gaussian import summarize_gaussian_log
from merlino_gic import (
    GICDefinition,
    define_gics_from_cartesian,
    evaluate_gic_definition,
    run_gicforge,
    write_gaussian_gic_input,
)
from merlino_gf import (
    run_gic_gf_report_from_fchk,
    run_gf_report_from_fchk,
    write_csv_tables as write_gf_csv_tables,
)
from merlino_semiexp import (
    DEFAULT_SEMIEXP_OBSERVABLE,
    DEFAULT_SEMIEXP_ROBUST_LOSS,
    DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
    HYDROGEN_PARAMETER_CONSTRAINT,
    ParameterClassConstraint,
    QMParameterPredicate,
    SemiexperimentalFitRequest,
    fit_semiexperimental_geometry,
    is_msr_legacy_file,
    read_observations,
    read_semiexperimental_job,
    semiexperimental_latex_tables,
    write_semiexperimental_html_report,
)
from merlino_vpt2_vci import (
    QuarticForceField,
    VCIOptions,
    load_force_field,
    run_vpt2_vci_report,
    solve_vci,
    write_csv_tables as write_vpt2_vci_csv_tables,
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

    gic_gf = sub.add_parser(
        "gic-gf",
        help="Run GF/PED from FCHK Hessian using a frozen GIC definition and optional Pulay scaling",
    )
    gic_gf.add_argument("--schema", type=Path, required=True, help="GIC definition JSON from gic-define")
    gic_gf.add_argument("--fchk", type=Path, required=True, help="Gaussian FCHK adapter containing Cartesian Hessian")
    gic_gf.add_argument("--geometry", "--xyz", dest="geometry", type=Path, help="Optional current Cartesian geometry for B")
    gic_gf.add_argument("--scale-file", type=Path, help="Optional Pulay diagonal scaling factors")
    gic_gf.add_argument(
        "--scale",
        action="append",
        default=[],
        help="Inline Pulay factor selector=value; selectors may be GIC001, name, one-based index, default or all",
    )
    gic_gf.add_argument("--out", type=Path)
    gic_gf.add_argument("--run-dir", type=Path)
    gic_gf.add_argument("--csv-dir", type=Path)

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
    gic.add_argument("--no-symmetry", action="store_true", help="Skip post-GICForge symmetry adaptation")

    gic_define = sub.add_parser("gic-define", help="Build a frozen GIC definition from Cartesian geometry")
    gic_define.add_argument("--geometry", "--xyz", dest="geometry", type=Path, required=True)
    gic_define.add_argument("--out", type=Path, required=True, help="Output GIC definition JSON")
    gic_define.add_argument("--workdir", type=Path, help="Optional GICForge working directory")
    gic_define.add_argument("--executable", type=Path)
    gic_define.add_argument("--gaussian-out", type=Path, help="Optional Gaussian-readable GIC block")
    gic_define.add_argument("--no-symmetry", action="store_true", help="Freeze raw non-redundant GICs without symmetry adaptation")

    gic_bmat = sub.add_parser("gic-bmatrix", help="Evaluate B matrix from a frozen GIC definition and Cartesian geometry")
    gic_bmat.add_argument("--schema", type=Path, required=True, help="GIC definition JSON from gic-define")
    gic_bmat.add_argument("--geometry", "--xyz", dest="geometry", type=Path, required=True)
    gic_bmat.add_argument("--out", type=Path, required=True, help="Output CSV B matrix")
    gic_bmat.add_argument("--values-out", type=Path, help="Optional CSV GIC values")
    gic_bmat.add_argument("--metadata-out", type=Path, help="Optional CSV GIC name/irrep metadata")

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
    semiexp.add_argument(
        "--job",
        type=Path,
        help="Merlino semiexperimental job file (.mfit, .mse.toml, .semiexp.toml) or legacy MSR file (.msr, .msr.inp)",
    )
    semiexp.add_argument(
        "--xyz",
        "--geometry",
        dest="xyz",
        type=Path,
        help="Initial parent Cartesian geometry in XYZ or Gaussian .com/.gjf format",
    )
    semiexp.add_argument(
        "--observations",
        type=Path,
        help="CSV/JSON/TOML with isotopologue B0 constants and corrections, or legacy MSR file (.msr, .msr.inp)",
    )
    semiexp.add_argument("--outdir", type=Path, required=True, help="Output directory for geometry, parameters, residuals and manifest")
    semiexp.add_argument("--backend", choices=("python", "fortran77"), default="python", help="Numerical backend requested by CLI/GUI")
    semiexp.add_argument(
        "--fixed",
        default="",
        help="Comma/semicolon-separated fixed GIC label patterns, primitive constraints or Gaussian-style GIC constraints such as NAME(Frozen,Value=0.0)=R[1,3]-R[1,2]",
    )
    semiexp.add_argument(
        "--fix-hydrogens",
        action="store_true",
        help="Freeze deterministic local H/D/T geometry constraints, expanded by symmetry",
    )
    semiexp.add_argument(
        "--max-iter",
        type=int,
        default=None,
        help="Maximum LM iterations; default is automatic: max(8, 2*N optimized parameters)",
    )
    semiexp.add_argument("--step", type=float, default=1.0e-4, help="Finite step for fallback derivatives with respect to working coordinates")
    semiexp.add_argument("--damping", type=float, default=1.0e-8, help="Initial Levenberg-Marquardt damping")
    semiexp.add_argument("--max-step", type=float, default=0.25, help="Maximum active-coordinate step norm per iteration")
    semiexp.add_argument(
        "--prune-condition",
        type=float,
        default=0.0,
        help="Auto-prune weak SE parameters until the initial weighted Jacobian condition is below this target; default 0 disables pruning",
    )
    semiexp.add_argument(
        "--robust-loss",
        choices=("none", "huber", "soft_l1", "cauchy"),
        default=DEFAULT_SEMIEXP_ROBUST_LOSS,
        help="Optional robust IRLS loss for experimental outlier isotopologues; default none",
    )
    semiexp.add_argument(
        "--robust-scale",
        type=float,
        default=0.0,
        help="Robust residual scale in weighted units; 0 selects automatic MAD scale",
    )
    semiexp.add_argument(
        "--leave-one-out",
        action="store_true",
        help="Run exact leave-one-isotopologue-out refits after the final SE fit",
    )
    semiexp.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help="Checkpoint path; default is semiexp_checkpoint.json in --outdir",
    )
    semiexp.add_argument(
        "--restart",
        type=Path,
        default=None,
        help="Restart from a Merlino SEfit checkpoint JSON",
    )
    semiexp.add_argument(
        "--observable",
        choices=("moments", "rotational_constants", "auto"),
        default=DEFAULT_SEMIEXP_OBSERVABLE,
        help="Fit target; default moments is the Merlino standard because it is more stable than reciprocal rotational constants",
    )
    semiexp.add_argument(
        "--coordinate-model",
        choices=("gic", "cartesian_symmetry"),
        default="gic",
        help="Working coordinates for SEfit: frozen GICs or symmetry-adapted Cartesians",
    )
    semiexp.add_argument(
        "--rotational-components",
        choices=("auto", "ABC", "AB", "AC", "BC"),
        default=DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
        help=(
            "Rotational component pair to use: planar molecules use only AB, AC or BC "
            "(mapped to moment pairs when observable=moments); auto chooses the most stable planar pair"
        ),
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

    semiexp_benchmark = sub.add_parser(
        "semiexp-benchmark",
        help="Generate versioned SEfit benchmark artifacts for regression tests and manuscripts",
    )
    semiexp_benchmark.add_argument("--paper", action="store_true", help="Generate the benchmark tables used in the paper")
    semiexp_benchmark.add_argument(
        "--snapshot",
        type=Path,
        default=Path("benchmarks/semiexp_msr/golden/semiexp_paper_regression.json"),
        help="Golden paper-regression snapshot",
    )
    semiexp_benchmark.add_argument(
        "--outdir",
        type=Path,
        default=Path("benchmarks/semiexp_msr/generated"),
        help="Directory for generated CSV and LaTeX benchmark tables",
    )
    semiexp_benchmark.add_argument(
        "--no-refresh",
        action="store_true",
        help="Use the snapshot exactly as written instead of refreshing available CSV outputs",
    )
    semiexp_benchmark.add_argument(
        "--update-snapshot",
        action="store_true",
        help="Overwrite the golden snapshot after refreshing from available CSV outputs",
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
            outputs.update({f"csv_{name}": path for name, path in write_gf_csv_tables(report, args.csv_dir).items()})
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

    if args.command == "gic-gf":
        report = run_gic_gf_report_from_fchk(
            args.fchk,
            args.schema,
            geometry_path=args.geometry,
            scale_path=args.scale_file,
            scale_records=tuple(args.scale),
        )
        out = args.out or Path("gic_gf_ped_report.txt")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report.text + "\n", encoding="utf-8")
        outputs = {"report": out}
        if args.csv_dir is not None:
            outputs.update({f"csv_{name}": path for name, path in write_gf_csv_tables(report, args.csv_dir, prefix="gic_gf").items()})
        run_dir = args.run_dir or out.parent
        inputs = {"fchk": args.fchk, "gic_definition": args.schema}
        if args.geometry is not None:
            inputs["geometry"] = args.geometry
        if args.scale_file is not None:
            inputs["scale_file"] = args.scale_file
        build_run_manifest(
            workflow="gic_gf",
            status="completed",
            run_dir=run_dir,
            inputs=inputs,
            outputs=outputs,
            parameters={"inline_scale_records": tuple(args.scale)},
            backend={"adapter": "gaussian-fchk", "solver": "python", "coordinate_model": "frozen-gic-definition"},
        ).write(Path(run_dir) / "gic_gf_manifest.json")
        print(out)
        print(f"frequency_count: {len(report.result.frequencies_cm)}")
        print(f"gic_count: {len(report.result.gic_labels)}")
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
            outputs.update({f"csv_{name}": path for name, path in write_vpt2_vci_csv_tables(report, args.csv_dir).items()})
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
        result = run_gicforge(args.workdir, executable=args.executable, symmetrize=not args.no_symmetry)
        print(f"manifest: {result.manifest}")
        print(f"symmetrized: {not args.no_symmetry}")
        for name, path in sorted(result.files.items()):
            print(f"{name}: {path}")
        return 0

    if args.command == "gic-define":
        from merlino_semiexp.geometry_input import read_geometry_input

        geometry = read_geometry_input(args.geometry)
        definition = define_gics_from_cartesian(
            tuple(geometry.atoms),
            geometry.coordinates_angstrom,
            workdir=args.workdir,
            executable=args.executable,
            symmetrize=not args.no_symmetry,
        )
        definition.write(args.out)
        print(f"schema: {args.out}")
        print(f"point_group: {definition.point_group}")
        print(f"symmetrized: {definition.symmetrized}")
        print(f"gic_count: {definition.u_matrix.shape[1]}")
        print(f"primitive_count: {definition.u_matrix.shape[0]}")
        if args.gaussian_out is not None:
            gaussian_path = write_gaussian_gic_input(definition, args.gaussian_out)
            print(f"gaussian_input: {gaussian_path}")
        return 0

    if args.command == "gic-bmatrix":
        import numpy as np
        from merlino_semiexp.geometry_input import read_geometry_input
        from topology.elements import atomic_number

        definition = GICDefinition.read(args.schema)
        geometry = read_geometry_input(args.geometry)
        z_numbers = tuple(atomic_number(atom) for atom in geometry.atoms)
        evaluation = evaluate_gic_definition(definition, geometry.coordinates_angstrom, atomic_numbers=z_numbers)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.out, evaluation.b_matrix, delimiter=",", fmt="%.16e")
        print(f"b_matrix: {args.out}")
        print(f"shape: {evaluation.b_matrix.shape[0]}x{evaluation.b_matrix.shape[1]}")
        print(f"point_group: {evaluation.point_group}")
        print(f"symmetrized: {evaluation.symmetrized}")
        if args.values_out is not None:
            args.values_out.parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(args.values_out, evaluation.values.reshape(-1, 1), delimiter=",", fmt="%.16e")
            print(f"values: {args.values_out}")
        if args.metadata_out is not None:
            args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
            rows = ["gic,name,irrep,label"]
            for idx, label in enumerate(evaluation.labels, start=1):
                name = evaluation.names[idx - 1] if idx <= len(evaluation.names) else f"GIC{idx:03d}"
                irrep = evaluation.irreps[idx - 1] if idx <= len(evaluation.irreps) else "UNK"
                rows.append(f"GIC{idx:03d},{name},{irrep},{json.dumps(label)}")
            args.metadata_out.write_text("\n".join(rows) + "\n", encoding="utf-8")
            print(f"metadata: {args.metadata_out}")
        return 0

    if args.command == "semiexp":
        legacy_msr_job = bool(args.job and is_msr_legacy_file(args.job))
        job = None if legacy_msr_job or not args.job else read_semiexperimental_job(args.job)
        geometry_path = args.xyz or (job.path if job is not None else None)
        observations_path = args.observations or (job.observations if job is not None else None)
        if legacy_msr_job:
            geometry_path = args.xyz or args.job
            observations_path = args.observations or args.job
        if geometry_path is None:
            raise ValueError("semiexp needs --geometry or --job")
        if observations_path is None:
            raise ValueError("semiexp needs --observations or a [files].observations entry in --job")
        fixed = _merge_unique(job.fixed_parameters if job else (), _parse_fixed_parameters(args.fixed))
        if args.fix_hydrogens:
            fixed = _merge_unique(fixed, (HYDROGEN_PARAMETER_CONSTRAINT,))
        observations = read_observations(observations_path)
        observable = _job_default(args.observable, DEFAULT_SEMIEXP_OBSERVABLE, job.observable if job else None)
        coordinate_model = _job_default(args.coordinate_model, "gic", job.coordinate_model if job else None)
        rotational_components = _job_default(
            args.rotational_components,
            DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
            job.rotational_components if job else None,
        )
        qm_predicates = _merge_unique(job.qm_predicates if job else (), _parse_qm_predicates(args.qm_predicate))
        parameter_classes = _merge_unique(job.parameter_classes if job else (), _parse_parameter_classes(args.parameter_class))
        backend = _job_default(args.backend, "python", job.backend if job else None)
        max_iter = args.max_iter if args.max_iter is not None else (job.max_iter if job else None)
        step = _job_default(args.step, 1.0e-4, job.step if job else None)
        damping = _job_default(args.damping, 1.0e-8, job.damping if job else None)
        max_step = _job_default(args.max_step, 0.25, job.max_step if job else None)
        prune_condition = _job_default(args.prune_condition, 0.0, job.prune_condition if job else None)
        robust_loss = _job_default(args.robust_loss, DEFAULT_SEMIEXP_ROBUST_LOSS, job.robust_loss if job else None)
        robust_scale = _job_default(args.robust_scale, 0.0, job.robust_scale if job else None)
        leave_one_out = bool(args.leave_one_out or (job.leave_one_out if job else False))
        checkpoint = args.checkpoint if args.checkpoint is not None else (job.checkpoint if job else None)
        restart = args.restart if args.restart is not None else (job.restart if job else None)
        request = SemiexperimentalFitRequest(
            initial_geometry=geometry_path,
            observations=observations,
            fixed_parameters=fixed,
            observable=observable,
            rotational_components=rotational_components,
            qm_predicates=qm_predicates,
            parameter_classes=parameter_classes,
            coordinate_model=coordinate_model,
            robust_loss=robust_loss,
            robust_scale=robust_scale,
            leave_one_out=leave_one_out,
        )
        result = fit_semiexperimental_geometry(
            request,
            max_iter=max_iter,
            step=step,
            damping=damping,
            max_step=max_step,
            prune_condition=prune_condition,
            checkpoint=checkpoint,
            restart=restart,
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
        rms_label = "rms_MHz" if result.diagnostics.observable == "rotational_constants" else "rms_observable"
        print(f"{rms_label}: {result.rms_MHz:.8g}")
        rot_diffs = [row.difference_MHz for row in result.rotational_constants]
        rotational_rms = math.sqrt(sum(diff * diff for diff in rot_diffs) / len(rot_diffs)) if rot_diffs else 0.0
        rotational_mse = sum(diff * diff for diff in rot_diffs) / len(rot_diffs) if rot_diffs else 0.0
        print(f"rotational_rms_MHz: {rotational_rms:.8g}")
        print(f"rotational_mean_square_MHz2: {rotational_mse:.8g}")
        print(f"rotational_mean_square_1e3_MHz2: {1000.0 * rotational_mse:.8g}")
        print(f"iterations: {result.iterations}")
        print(f"stationary_point: {result.stationary_point}")
        print(f"convergence: {result.diagnostics.convergence_reason}")
        print(f"rank: {result.diagnostics.rank}")
        print(f"condition_number: {result.diagnostics.condition_number:.8g}")
        print(f"observable: {result.diagnostics.observable}")
        print(f"components: {','.join(result.diagnostics.components)}")
        print(f"backend: {backend}")
        print(f"coordinate_model: {result.diagnostics.coordinate_model}")
        return 0

    if args.command == "semiexp-benchmark":
        if not args.paper:
            raise ValueError("semiexp-benchmark currently requires --paper")
        from merlino_semiexp.paper_benchmarks import generate_paper_benchmark_artifacts

        _snapshot, artifacts = generate_paper_benchmark_artifacts(
            snapshot_path=args.snapshot,
            outdir=args.outdir,
            refresh_from_outputs=not args.no_refresh,
            update_snapshot=args.update_snapshot,
        )
        for name, path in artifacts.items():
            print(f"{name}: {path}")
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
    return tuple(part.strip() for part in _split_top_level(raw, separators=",;") if part.strip())


def _split_top_level(raw: str, *, separators: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    round_depth = 0
    square_depth = 0
    brace_depth = 0
    for char in str(raw):
        if char == "(":
            round_depth += 1
        elif char == ")" and round_depth > 0:
            round_depth -= 1
        elif char == "[":
            square_depth += 1
        elif char == "]" and square_depth > 0:
            square_depth -= 1
        elif char == "{":
            brace_depth += 1
        elif char == "}" and brace_depth > 0:
            brace_depth -= 1
        if char in separators and round_depth == 0 and square_depth == 0 and brace_depth == 0:
            parts.append("".join(current))
            current = []
            continue
        current.append(char)
    parts.append("".join(current))
    return parts


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


def _merge_unique[T](left: tuple[T, ...], right: tuple[T, ...]) -> tuple[T, ...]:
    result: list[T] = []
    for item in (*left, *right):
        if item not in result:
            result.append(item)
    return tuple(result)


def _job_default[T](value: T, default: T, job_value: T | None) -> T:
    if job_value is not None and value == default:
        return job_value
    return value


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
