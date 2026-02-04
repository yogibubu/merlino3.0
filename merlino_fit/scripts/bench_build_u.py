#!/usr/bin/env python
from __future__ import annotations

import argparse
import time
import numpy as np

from survibfit.modify_geom import read_xyz
from survibfit.pipeline import primitives_from_topology, build_topology
from survibfit.transforms import build_u
from survibfit.pipeline import _load_topology_elements


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--repeat", type=int, default=5)
    ap.add_argument("--symmetry-mode", default="hybrid")
    ap.add_argument("--prune-mode", default="svd")
    ap.add_argument("--zeff-tol", type=float, default=0.05)
    ap.add_argument("--geometry-match-tol", type=float, default=12.0)
    args = ap.parse_args()

    atoms, coords_ang, _ = read_xyz(args.xyz)
    coords_au = coords_ang / 0.52917721092
    atomic_number = _load_topology_elements()
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)

    prims = primitives_from_topology(coords_au, Z, np.deg2rad(170.0))
    _, _, ringset = build_topology(coords_au, Z)

    t0 = time.perf_counter()
    for _ in range(args.repeat):
        build_u(
            prims,
            coords_au,
            Z=Z,
            ringset=ringset,
            symmetry_mode=args.symmetry_mode,
            prune_mode=args.prune_mode,
            zeff_tol=args.zeff_tol,
            geometry_match_tol=args.geometry_match_tol,
        )
    t1 = time.perf_counter()
    print(f"build_u avg: {(t1 - t0) / args.repeat:.6f} s")


if __name__ == "__main__":
    main()
