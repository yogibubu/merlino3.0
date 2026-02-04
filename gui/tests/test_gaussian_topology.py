from pathlib import Path

from gui.gaussian import _build_gaussian_topology_lines
from topology.test_topology import _parse_gaussian_topology_overrides


def test_build_gaussian_topology_lines_parses_cm5_and_mayer():
    lines = [
        " Hirshfeld charges, spin densities, dipoles, and CM5 charges",
        "   1  C   -0.0100   0.0000  -0.1234",
        "   2  H    0.0100   0.0000   0.0456",
        "",
        " Mayer bond orders and valences:",
        "  B( 1-C, 2-H)=0.9123",
    ]
    out = _build_gaussian_topology_lines(lines)
    assert out is not None
    joined = "\n".join(out)
    assert "CM5 1" in joined
    assert "CM5 2" in joined
    assert "BO_SOURCE = Mayer" in joined
    assert "BO 1 2" in joined


def test_build_gaussian_topology_lines_parses_mayer_atomic_matrix():
    lines = [
        " Atomic Valencies and Mayer Atomic Bond Orders:",
        "              1          2          3",
        "    1  C    3.620819   1.310133  -0.075588",
        "    2  C    1.310133   3.202272   1.197745",
        "    3  C   -0.075588   1.197745   3.569945",
    ]
    out = _build_gaussian_topology_lines(lines)
    assert out is not None
    joined = "\n".join(out)
    assert "BO_SOURCE = Mayer" in joined
    assert "BO 1 2" in joined
    assert "BO 2 3" in joined
    assert "BO 1 3" not in joined


def test_parse_gaussian_topology_overrides_section(tmp_path: Path):
    xyzin = tmp_path / "xyzin"
    xyzin.write_text(
        "2\n"
        "x\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 0.7\n"
        "#GAUSSIAN_TOPOLOGY\n"
        "CM5_COUNT = 2\n"
        "CM5 1 -0.1000\n"
        "CM5 2  0.1000\n"
        "BO_SOURCE = Mayer\n"
        "BO_COUNT = 1\n"
        "BO 1 2 0.8500\n"
    )
    cm5, bo, meta = _parse_gaussian_topology_overrides(xyzin)
    assert cm5 == {0: -0.1, 1: 0.1}
    assert bo == {(0, 1): 0.85}
    assert meta["charge_source"] == "Gaussian CM5"
    assert "Gaussian" in meta["bond_order_source"]
