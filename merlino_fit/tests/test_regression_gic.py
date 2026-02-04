from pathlib import Path
import numpy as np
import re

from survibfit.modify_geom import read_xyz
from merlino_fit.topology.elements import atomic_number
from survibfit.pipeline import build_topology, primitives_from_topology
from survibfit.transforms import build_u_with_names, format_readgic_lines


def _generate_lines(xyz_path: Path):
    atoms, coords_ang, _ = read_xyz(xyz_path)
    coords_au = coords_ang / 0.52917721092
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    prims = primitives_from_topology(coords_au, Z, np.deg2rad(170.0))
    _, _, ringset = build_topology(coords_au, Z)
    _, names = build_u_with_names(
        prims,
        coords_au,
        Z=Z,
        ringset=ringset,
        include_frag=False,
        symmetry_mode="gblock",
    )
    return format_readgic_lines(names)


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
