import numpy as np

from survibfit.pipeline import primitives_from_topology
from survibfit.primitives import Primitive


def test_fragment_primitives_stable():
    # Two fragments (H2O + H2O) in angstrom
    coords = np.array(
        [
            [0.0, 0.0, 0.0],   # O
            [0.96, 0.0, 0.0],  # H
            [-0.24, 0.93, 0.0],  # H
            [0.0, 0.0, 3.0],   # O
            [0.96, 0.0, 3.0],  # H
            [-0.24, 0.93, 3.0],  # H
        ],
        dtype=float,
    )
    Z = np.array([8, 1, 1, 8, 1, 1], dtype=int)
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    frag = [p for p in prims if p.kind.startswith("frag_")]
    assert len(frag) == 6


def test_dihedral_u_pick_and_combine():
    # Branching around bond (1-2) to create multiple dihedrals
    coords = np.array(
        [
            [0.0, 0.0, 0.0],    # 0 (C)
            [1.54, 0.0, 0.0],   # 1 (C)
            [-0.6, 1.0, 0.2],   # 2 (C)
            [-0.6, -1.0, -0.2], # 3 (C)
            [2.1, 1.0, -0.2],   # 4 (C)
            [2.1, -1.0, 0.2],   # 5 (C)
        ],
        dtype=float,
    )
    Z = np.array([6, 6, 6, 6, 6, 6], dtype=int)
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    from survibfit.transforms import dihedral_u

    U_pick, idx_pick = dihedral_u(prims, coords, Z=Z, priority="zeff", mode="pick", coords_units="ang")
    U_comb, idx_comb = dihedral_u(prims, coords, Z=Z, priority="zeff", mode="combine", coords_units="ang")

    assert U_pick.shape[0] == len(idx_pick)
    assert U_comb.shape[0] == len(idx_comb)
    # combine should reduce each multi-dihedral group to a single column
    assert U_comb.shape[1] <= U_pick.shape[1]


def test_dihedral_u_orbit_closure(monkeypatch):
    from survibfit import transforms
    from survibfit.transforms import dihedral_u

    prims = [
        Primitive("bond", (0, 1)),
        Primitive("bond", (4, 1)),
        Primitive("bond", (1, 2)),
        Primitive("bond", (2, 3)),
        Primitive("bond", (2, 5)),
        Primitive("dihedral", (0, 1, 2, 3)),
        Primitive("dihedral", (0, 1, 2, 5)),
        Primitive("dihedral", (4, 1, 2, 3)),
        Primitive("dihedral", (4, 1, 2, 5)),
    ]
    coords = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [2.5, 0.5, 0.2],
            [-0.8, 1.0, 0.1],
            [2.4, -0.4, 0.3],
        ],
        dtype=float,
    )
    Z = np.array([6, 6, 6, 6, 6, 6], dtype=int)
    priority = np.array([1.0, 2.0, 2.0, 1.5, 3.0, 3.0], dtype=float)

    def fake_atom_classes(coords_in, Z_in, zeff_tol=0.05):
        atom_class = np.array([10, 20, 30, 40, 10, 40], dtype=int)
        return atom_class, {}, None, None, None

    monkeypatch.setattr(transforms, "atom_classes", fake_atom_classes)

    U, idx, info = dihedral_u(
        prims,
        coords,
        Z=Z,
        priority=priority,
        mode="orbit",
        return_info=True,
    )

    assert idx == [5, 6, 7, 8]
    assert U.shape == (4, 1)
    assert np.allclose(U[:, 0], np.full(4, 0.25))
    assert info["bonds"][0]["reason"] == "orbit"
    assert info["bonds"][0]["selected"] == [5, 6, 7, 8]


def test_xy3_valence_block_matches_fortran_style():
    from survibfit.transforms import valence_angle_u

    coords = np.array(
        [
            [0.0, 0.0, 0.0],   # center
            [1.0, 0.0, 0.0],   # heavy
            [-0.5, 0.87, 0.0], # heavy
            [0.0, -1.0, 0.0],  # H
        ],
        dtype=float,
    )
    Z = np.array([7, 6, 6, 1], dtype=int)
    prims = [
        Primitive("bond", (0, 1)),
        Primitive("bond", (0, 2)),
        Primitive("bond", (0, 3)),
        Primitive("angle", (1, 0, 2)),
        Primitive("angle", (1, 0, 3)),
        Primitive("angle", (2, 0, 3)),
    ]
    U, idx = valence_angle_u(prims, coords, ringset=None, Z=Z)

    angle_map = {tuple(sorted((prims[i].atoms[0], prims[i].atoms[2]))): r for r, i in enumerate(idx)}
    # Expected angle primitive order from the explicit XY3 basis:
    # (1,0,2), (1,0,3), (2,0,3)
    assert U.shape[1] == 2
    assert np.allclose(U[angle_map[(1, 2)], :], [2.0 / np.sqrt(6.0), 0.0], atol=1e-8)
    assert np.allclose(U[angle_map[(1, 3)], :], [-1.0 / np.sqrt(6.0), 1.0 / np.sqrt(2.0)], atol=1e-8)
    assert np.allclose(U[angle_map[(2, 3)], :], [-1.0 / np.sqrt(6.0), -1.0 / np.sqrt(2.0)], atol=1e-8)


def test_xy2_valence_block_is_explicit():
    from survibfit.transforms import valence_angle_u

    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    Z = np.array([8, 1, 1], dtype=int)
    prims = [
        Primitive("bond", (0, 1)),
        Primitive("bond", (0, 2)),
        Primitive("angle", (1, 0, 2)),
    ]
    U, idx = valence_angle_u(prims, coords, ringset=None, Z=Z)
    assert idx == [2]
    assert U.shape == (1, 1)
    assert np.allclose(U, np.eye(1))


def test_ring_angle_u_excludes_valence():
    # Simple 4-membered ring (square) to trigger ring angles
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [1.5, 1.5, 0.0],
            [0.0, 1.5, 0.0],
        ],
        dtype=float,
    )
    Z = np.array([6, 6, 6, 6], dtype=int)
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    from survibfit.pipeline import build_topology
    from survibfit.transforms import valence_angle_u, ring_angle_u

    _, _, ringset = build_topology(coords, Z)
    U_val, idx_val = valence_angle_u(prims, coords, ringset=ringset)
    U_ring, idx_ring = ring_angle_u(prims, coords, ringset=ringset)

    # All angles in the 4-ring should be handled by ring U, not valence U
    assert len(idx_ring) > 0
    assert len(idx_val) == 0


def test_oop_u_skips_near_linear():
    # Three neighbors nearly linear around center should drop OOPs
    coords = np.array(
        [
            [0.0, 0.0, 0.0],     # center
            [1.0, 0.0, 0.0],     # neighbor
            [-1.0, 0.0, 0.0],    # neighbor (linear)
            [0.0, 1.0, 0.0],     # neighbor
        ],
        dtype=float,
    )
    Z = np.array([6, 8, 8, 8], dtype=int)
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    from survibfit.transforms import oop_u

    U_oop, idx_oop = oop_u(prims, coords, linear_threshold=np.deg2rad(170.0))
    # near-linear center (atom 0) => no oop kept with center 0
    for idx in idx_oop:
        _, j, _, _ = prims[idx].atoms
        assert j != 0


def test_build_u_dimensions():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [1.8, 1.0, 0.0],
            [2.8, 1.0, 0.5],
        ],
        dtype=float,
    )
    Z = np.array([6, 6, 6, 6], dtype=int)
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    from survibfit.transforms import build_u
    from survibfit.pipeline import build_topology

    _, _, ringset = build_topology(coords, Z)
    U = build_u(prims, coords, Z=Z, ringset=ringset)
    assert U.shape[0] == len(prims)
    assert U.shape[1] > 0


def test_build_u_with_names():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    Z = np.array([8, 8, 8], dtype=int)
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    from survibfit.transforms import build_u_with_names, format_readgic_lines
    from survibfit.pipeline import build_topology

    _, _, ringset = build_topology(coords, Z)
    U, names = build_u_with_names(prims, coords, Z=Z, ringset=ringset, include_frag=False)
    lines = format_readgic_lines(names)
    assert U.shape[1] == len(names)
    assert all("(" in line and ")" in line for line in lines)
