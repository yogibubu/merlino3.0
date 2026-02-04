from pathlib import Path

import numpy as np

from topology.pipeline import build_topology_objects
from topology.topology_reporting import print_topology_report


def test_topology_report_contains_symmetry_and_equivalent_classes(tmp_path, monkeypatch):
    # Idealized methane geometry (angstrom)
    a = 1.089
    u = a / np.sqrt(3.0)
    coords = np.array(
        [
            [0.0, 0.0, 0.0],      # C
            [u, u, u],            # H
            [u, -u, -u],          # H
            [-u, u, -u],          # H
            [-u, -u, u],          # H
        ],
        dtype=float,
    )
    Z = [6, 1, 1, 1, 1]

    cg, dg, _ringset, synthons, aromaticity = build_topology_objects(coords, Z)

    monkeypatch.chdir(tmp_path)
    print_topology_report(cg, dg, synthons, arom=aromaticity, filename="topology.report")

    report_path = Path(tmp_path) / "working" / "topology.report"
    text = report_path.read_text()

    assert "GLOBAL SYMMETRY" in text
    assert "Point group:" in text
    assert any(pg in text for pg in ("Td", "T", "Th", "C3v"))
    assert "EQUIVALENT INTERNAL-PARAMETER CLASSES" in text
    assert "bond:" in text
    assert "Equivalent atom classes" in text
