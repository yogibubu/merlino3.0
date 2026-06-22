from pathlib import Path
import re
import subprocess
import sys
import tempfile

import pytest

from survibfit.modify_geom import read_xyz
from merlino_gic import define_gics_from_cartesian
from merlino_gic.model import parse_gicforge_line


def _generate_lines(xyz_path: Path):
    atoms, coords_ang, _ = read_xyz(xyz_path)
    with tempfile.TemporaryDirectory(prefix="survibfit_gic_regression_") as tmp:
        definition = define_gics_from_cartesian(tuple(atoms), coords_ang, workdir=Path(tmp), symmetrize=False)
        return [line.rstrip() for line in definition.gaussian_input.splitlines() if parse_gicforge_line(line) is not None]


def _flip_linear_combo_sign(line: str) -> str:
    """Flip global sign for Value and all coefficients in bracket expressions."""
    if "=[" not in line:
        return line

    out = line
    m = re.search(r"Value=\s*([+-]?\d+\.\d+)", out)
    if m:
        val = float(m.group(1))
        out = out[: m.start(1)] + f"{-val: .5f}" + out[m.end(1) :]

    def _flip_coeff(mm):
        c = float(mm.group(1))
        return f"{-c: .5f}*"

    out = re.sub(r"([+-]?\d+\.\d+)\*", _flip_coeff, out)
    return out


def _assert_lines_sign_robust(expected, got):
    def norm(s: str) -> str:
        s = " ".join(s.split())
        def _norm_val(mm):
            return f"Value={float(mm.group(1)):.5f}"
        s = re.sub(r"Value=\s*([+-]?\d+\.\d+)", _norm_val, s)
        return s

    assert len(expected) == len(got)
    for e, g in zip(expected, got):
        if norm(e) == norm(g):
            continue
        assert norm(e) == norm(_flip_linear_combo_sign(g))


def test_gic_regression_ooo_triangle():
    base = Path(__file__).resolve().parent
    xyz = base / "data" / "ooo_triangle.xyz"
    golden = base / "golden" / "ooo_triangle.gic.txt"
    lines = _generate_lines(xyz)
    _assert_lines_sign_robust(
        [l for l in golden.read_text().splitlines() if l.strip()],
        lines,
    )


def test_gic_regression_c4_chain():
    base = Path(__file__).resolve().parent
    xyz = base / "data" / "c4_chain.xyz"
    golden = base / "golden" / "c4_chain.gic.txt"
    lines = _generate_lines(xyz)
    _assert_lines_sign_robust(
        [l for l in golden.read_text().splitlines() if l.strip()],
        lines,
    )


def test_survibfit_gic_cli_is_gicforge_identical(tmp_path):
    base = Path(__file__).resolve().parent
    xyz = base / "data" / "c4_chain.xyz"
    out = tmp_path / "python.gic"
    workdir = tmp_path / "gicforge"

    subprocess.run(
        [
            sys.executable,
            "-m",
            "survibfit.cli",
            "gic",
            "--xyz",
            str(xyz),
            "--out",
            str(out),
            "--workdir",
            str(workdir),
        ],
        check=True,
    )

    python_lines = out.read_text(encoding="utf-8").splitlines()
    fortran_lines = [
        line.rstrip()
        for line in (workdir / "gauin").read_text(encoding="utf-8", errors="replace").splitlines()
        if parse_gicforge_line(line) is not None
    ]
    assert python_lines == fortran_lines


@pytest.mark.parametrize(
    "xyz_name",
    [
        "polycyclics/naphthalene.xyz",
        "polycyclics/anthracene.xyz",
        "polycyclics/coronene.xyz",
        "polycyclics/norbornane.xyz",
        "polycyclics/myrtenol.xyz",
        "polycyclics/testosterone.xyz",
        "polycyclics/saccharine.xyz",
        "polycyclics/2_deoxyribose.xyz",
    ],
)
def test_survibfit_gic_cli_is_gicforge_identical_polycyclic(tmp_path, xyz_name):
    base = Path(__file__).resolve().parent
    xyz = base / "data" / xyz_name
    out = tmp_path / f"{Path(xyz_name).stem}.gic"
    workdir = tmp_path / Path(xyz_name).stem

    subprocess.run(
        [
            sys.executable,
            "-m",
            "survibfit.cli",
            "gic",
            "--xyz",
            str(xyz),
            "--out",
            str(out),
            "--workdir",
            str(workdir),
        ],
        check=True,
    )

    python_lines = out.read_text(encoding="utf-8").splitlines()
    fortran_lines = [
        line.rstrip()
        for line in (workdir / "gauin").read_text(encoding="utf-8", errors="replace").splitlines()
        if parse_gicforge_line(line) is not None
    ]
    assert python_lines == fortran_lines
