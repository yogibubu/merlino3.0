import numpy as np

from survibfit.pipeline import primitives_from_topology, b_matrix, build_topology
from survibfit.transforms import valence_angle_u, oop_u


def test_valence_angle_u_rank():
    # center + three neighbors
    coords = np.array(
        [
            [0.0, 0.0, 0.0],   # center
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    Z = np.array([6, 8, 8, 8], dtype=int)
    coords_au = coords / 0.52917721092
    prims = primitives_from_topology(coords_au, Z, linear_threshold=np.deg2rad(170.0))
    _, _, ringset = build_topology(coords_au, Z)

    U, idx = valence_angle_u(prims, coords, ringset=ringset)
    assert len(idx) > 0

    prim_block = [prims[i] for i in idx]
    Bc = b_matrix(prim_block, coords_au, 1e-4)
    G = Bc @ Bc.T
    evals, _ = np.linalg.eigh(G)
    rank = int(np.sum(evals > 1e-8))
    assert U.shape[1] == rank
    # Check diagonalization
    M = U.T @ G @ U
    off = M - np.diag(np.diag(M))
    assert np.linalg.norm(off) < 1e-6


def test_oop_u_rank():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],   # center
            [1.0, 0.0, 0.0],
            [-0.5, 0.8, 0.0],
            [-0.5, -0.8, 0.0],
        ],
        dtype=float,
    )
    Z = np.array([6, 8, 8, 8], dtype=int)
    coords_au = coords / 0.52917721092
    prims = primitives_from_topology(coords_au, Z, linear_threshold=np.deg2rad(170.0))
    U, idx = oop_u(prims, coords_au)
    assert U.shape[1] <= len(idx)
