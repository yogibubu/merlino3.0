#!/usr/bin/env python3
"""
Integration test / pipeline driver: xyzin -> topology -> xyzin (append)
"""

import sys
from pathlib import Path
import numpy as np

from .elements import atomic_number
from .continuous_graph import bond_order
from .discrete_graph import DiscreteGraph
from .ringset import RingSet
from .atomic_synthons import AtomicSynthons
from .aromaticity import Aromaticity
from .topology_writer import write_topology_section


class ContinuousGraphAdapter:
    def __init__(self, coords, Z, bond_order_overrides=None):
        self.coords = coords
        self.Z = Z
        self.natoms = len(Z)
        self.bond_order_overrides = bond_order_overrides or {}

        neighbors = [list(range(self.natoms)) for _ in range(self.natoms)]
        for i in range(self.natoms):
            neighbors[i].remove(i)

        self.BO = np.zeros((self.natoms, self.natoms))
        cache = {}

        for i in range(self.natoms):
            for j in range(i + 1, self.natoms):
                bo = self.bond_order_overrides.get((i, j))
                if bo is None:
                    bo = bond_order(i, j, Z, coords, neighbors, cache)
                self.BO[i, j] = self.BO[j, i] = bo


def _parse_gaussian_topology_overrides(xyzin: Path):
    lines = xyzin.read_text().splitlines()
    try:
        i0 = next(i for i, ln in enumerate(lines) if ln.strip().upper() == "#GAUSSIAN_TOPOLOGY")
    except StopIteration:
        return {}, {}, {
            "charge_source": "Synthons model",
            "bond_order_source": "Topology continuous model",
        }

    cm5 = {}
    bo = {}
    bo_source = None
    for k in range(i0 + 1, len(lines)):
        s = lines[k].strip()
        if not s:
            continue
        if s.startswith("#"):
            break
        parts = s.split()
        if len(parts) >= 3 and parts[0].upper() == "CM5":
            try:
                idx = int(parts[1]) - 1
                val = float(parts[2])
                if idx >= 0:
                    cm5[idx] = val
            except Exception:
                continue
        elif len(parts) >= 3 and parts[0].upper() == "BO_SOURCE":
            bo_source = parts[2] if parts[1] == "=" else parts[1]
        elif len(parts) >= 4 and parts[0].upper() == "BO":
            try:
                i = int(parts[1]) - 1
                j = int(parts[2]) - 1
                val = float(parts[3])
                if i >= 0 and j >= 0 and i != j:
                    key = (i, j) if i < j else (j, i)
                    bo[key] = val
            except Exception:
                continue
    meta = {
        "charge_source": "Gaussian CM5" if cm5 else "Synthons model",
        "bond_order_source": (
            f"Gaussian {bo_source}" if bo and bo_source else
            ("Gaussian override" if bo else "Topology continuous model")
        ),
    }
    return cm5, bo, meta


def read_xyz_from_xyzin(xyzin):
    with open(xyzin) as fh:
        lines = fh.readlines()

    nat = int(lines[0])
    xyz = lines[2:2 + nat]

    atoms, coords = [], []
    for l in xyz:
        f = l.split()
        atoms.append(f[0])
        coords.append([float(x) for x in f[1:4]])

    Z = [atomic_number(a) for a in atoms]
    return np.array(coords), Z


def run_topology_on_xyzin(
    xyzin,
    *,
    symm_tol: float = 1.0e-3,
    symmetrize_coords: bool = False,
):
    from geometry.thermo_trasl import parse_xyzin_basic_section
    from .topology_reporting import print_topology_report

    try:
        from .rdkit_bridge import xyz_to_smiles
    except Exception:
        xyz_to_smiles = None

    xyzin = Path(xyzin)

    coords, Z = read_xyz_from_xyzin(xyzin)
    cm5_overrides, bo_overrides, meta = _parse_gaussian_topology_overrides(xyzin)

    cg = ContinuousGraphAdapter(coords, Z, bond_order_overrides=bo_overrides)
    dg = DiscreteGraph(cg)

    ringset = RingSet(dg, coords=cg.coords)

    neighbors = [list(dg.adjacency[i]) for i in range(dg.natoms)]

    synthons = AtomicSynthons(
        Z=Z,
        coords=cg.coords,
        neighbors=neighbors,
    )
    synthons._external_charges = cm5_overrides or None
    synthons._external_bond_orders = bo_overrides or None
    synthons._charge_source = meta["charge_source"]
    synthons._bond_order_source = meta["bond_order_source"]

    aromaticity = Aromaticity(
        graph=cg,
        discrete_graph=dg,
        ring_set=ringset,
        force_aromatic=False,
    )

    smiles = ""
    if xyz_to_smiles is not None:
        smiles = xyz_to_smiles(
            coords,
            Z,
            aromatic_smiles=True,
            force_aromatic=False,
        )

    rep = parse_xyzin_basic_section(str(xyzin)).get("REPRESENTATION", "Ir")
    print_topology_report(
        cg=cg,
        dg=dg,
        synthons=synthons,
        arom=aromaticity,
        ringset=ringset,
        filename="topology.report",
        symm_tol=float(symm_tol),
        symmetrize_coords=bool(symmetrize_coords),
    )

    with open(xyzin, "a") as fh:
        write_topology_section(
            fh,
            cg=cg,
            dg=dg,
            ringset=ringset,
            aromaticity=aromaticity,
            synthons=synthons,
            spin_density=None,
            smiles=smiles,
        )


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: test_topology.py xyzin")
    run_topology_on_xyzin(sys.argv[1])


if __name__ == "__main__":
    main()
