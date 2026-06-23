import numpy as np
from pathlib import Path

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


def test_real_polycyclic_cases_build_consistently():
    from survibfit.modify_geom import read_xyz
    from topology.elements import atomic_number
    from survibfit.transforms import build_u, build_u_blocks

    cases = {
        "naphthalene": ("data/polycyclics/naphthalene.xyz", 18, 167, 48),
        "anthracene": ("data/polycyclics/anthracene.xyz", 24, 238, 66),
        "phenanthrene": ("data/polycyclics/phenanthrene.xyz", 24, 238, 66),
        "pyrene": ("data/polycyclics/pyrene.xyz", 26, 277, 72),
        "coronene": ("data/polycyclics/coronene.xyz", 36, 426, 102),
        "dihydronaphthalene": ("data/polycyclics/dihydronaphthalene.xyz", 20, 209, 54),
        "norbornane": ("data/polycyclics/norbornane.xyz", 19, 290, 51),
        "norbornene": ("data/polycyclics/norbornene.xyz", 17, 242, 45),
        "norbornadiene": ("data/polycyclics/norbornadiene.xyz", 15, 194, 39),
        "myrtenol": ("data/polycyclics/myrtenol.xyz", 27, 361, 75),
        "testosterone": ("data/polycyclics/testosterone.xyz", 49, 723, 141),
        "saccharine": ("data/polycyclics/saccharine.xyz", 17, 172, 45),
        "2_deoxyribose": ("data/polycyclics/2_deoxyribose.xyz", 19, 215, 51),
    }
    for name, (relpath, nat_expected, nprim_expected, nactive_expected) in cases.items():
        path = Path(__file__).resolve().parent / relpath
        atoms, coords_ang, _ = read_xyz(path)
        Z = np.array([atomic_number(a) for a in atoms], dtype=int)
        coords = coords_ang / 0.52917721092
        prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))
        _, _, ringset = build_topology(coords, Z)
        U = build_u(prims, coords, Z=Z, ringset=ringset, symmetry_mode="hybrid", prune_mode="none")
        assert len(atoms) == nat_expected, name
        assert len(prims) == nprim_expected, name
        assert U.shape == (nprim_expected, nactive_expected), name
        assert U.shape[0] == len(prims)
        assert U.shape[1] > 0
        labels = [label for label, _, _ in build_u_blocks(prims, coords, Z=Z, ringset=ringset, symmetry_mode="hybrid", prune_mode="none")]
        assert "bond" in labels
        assert "dihedral" in labels
        assert any(label in labels for label in ("cyclic_valence_bend", "butterfly", "hinge"))
        assert "out_of_plane" in labels
