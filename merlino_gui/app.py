from __future__ import annotations

import argparse
import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication

from merlino_core import ensure_project_state, repo_root
from .dashboard import DashboardWindow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Launch the Merlino4 workflow dashboard.")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help="Project working directory. Defaults to <repo>/working.",
    )
    return parser


def run(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = repo_root(Path(__file__))
    workdir = args.workdir.expanduser() if args.workdir is not None else root / "working"
    state = ensure_project_state(workdir)

    app = QApplication(sys.argv[:1])
    window = DashboardWindow(state.workdir)
    window.show()
    return int(app.exec())


def main() -> None:
    raise SystemExit(run(sys.argv[1:]))


if __name__ == "__main__":
    main()
