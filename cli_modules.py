"""
cli.py
======

Command-line interface for molecular pipelines.

No pipeline is executed unless explicitly requested.
"""

import argparse
from pathlib import Path
import sys

from gui.project_manager import ProjectManager


# ============================================================
# CLI
# ============================================================
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Molecular pipeline driver"
    )

    parser.add_argument(
        "working_dir",
        type=Path,
        help="Working directory containing xyzin",
    )

    parser.add_argument(
        "--rotational",
        action="store_true",
        help="Run rotational (spectroscopic) pipeline",
    )

    parser.add_argument(
        "--topology",
        action="store_true",
        help="Run topology pipeline",
    )

    parser.add_argument(
        "--thermo",
        action="store_true",
        help="Run thermodynamic pipeline",
    )

    parser.add_argument(
        "--rovib",
        action="store_true",
        help="Run rovibrational DOS/Q(T) pipeline",
    )

    parser.add_argument("--rovib-vib-dos", help="Vibrational DOS file (default: dos_vib.dat)")
    parser.add_argument("--rovib-out", help="Rovibrational DOS output (default: dos_rovib.dat)")
    parser.add_argument("--rovib-rot-out", help="Rotational DOS output (optional)")
    parser.add_argument("--rovib-qout", help="Rovibrational Q(T) output (default: rovib_qt.dat)")
    parser.add_argument("--rovib-t", type=float, help="Temperature in K (default: T_K in #BASIC)")
    parser.add_argument("--rovib-emax-rot", type=float, help="Max rotational energy (cm-1)")
    parser.add_argument("--rovib-jmax", type=int, help="Max J for rotational levels")
    parser.add_argument("--rovib-sigma", type=int, help="Symmetry number override")
    parser.add_argument("--rovib-rotor-type", help="Rotor type override")
    parser.add_argument("--rovib-A", type=float, help="A rotational constant (MHz)")
    parser.add_argument("--rovib-B", type=float, help="B rotational constant (MHz)")
    parser.add_argument("--rovib-C", type=float, help="C rotational constant (MHz)")

    parser.add_argument(
        "--input-only",
        action="store_true",
        help="Validate input only, run no pipelines",
    )

    args = parser.parse_args(argv)

    pm = ProjectManager(args.working_dir)

    # --------------------------------------------------------
    # Always validate input (NO try/except: explode on error)
    # --------------------------------------------------------
    pm.validate_xyzin()

    # --------------------------------------------------------
    # Decide pipelines
    # --------------------------------------------------------
    if args.input_only:
        print("Input validated. No pipelines executed.")
        return 0

    pipelines = []
    if args.rotational:
        pipelines.append("rotational")
    if args.topology:
        pipelines.append("topology")
    if args.thermo:
        pipelines.append("thermo")
    if args.rovib:
        pipelines.append("rovib")

    if not pipelines:
        print(
            "Nothing to do. "
            "Specify at least one pipeline or --input-only.",
            file=sys.stderr,
        )
        return 1

    # --------------------------------------------------------
    # Run pipelines (NO try/except: explode on error)
    # --------------------------------------------------------
    pm.rovib_options = {
        "vib_dos": args.rovib_vib_dos,
        "out": args.rovib_out,
        "rot_out": args.rovib_rot_out,
        "q_out": args.rovib_qout,
        "t_k": args.rovib_t,
        "emax_rot": args.rovib_emax_rot,
        "jmax": args.rovib_jmax,
        "sigma": args.rovib_sigma,
        "rotor_type": args.rovib_rotor_type,
        "A_MHz": args.rovib_A,
        "B_MHz": args.rovib_B,
        "C_MHz": args.rovib_C,
    }

    pm.run(pipelines)

    print("Completed pipelines:", ", ".join(pipelines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
