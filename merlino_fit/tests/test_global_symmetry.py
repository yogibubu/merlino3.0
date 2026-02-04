import numpy as np
from pathlib import Path

from survibfit.modify_geom import read_xyz
from topology.elements import atomic_number
from survibfit.pipeline import primitives_from_topology
from survibfit.symmetry_global import orient_coords, symmetry_elements_from_geometry, symmetrize_u


def test_c5_axis_present():
    xyz = Path(__file__).resolve().parent / "data" / "c5h5.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    symbols = [a for a in atoms]
    coords = coords_ang

    coords_oriented = orient_coords(coords, weights=Z)
    elements, classes, permutations = symmetry_elements_from_geometry(
        symbols, coords_oriented, tol=1.0e-2, max_n=8
    )
    labels = [e[0] for e in elements]
    assert any(lab.startswith("C5") for lab in labels)


def test_a1_labels_present():
    xyz = Path(__file__).resolve().parent / "data" / "c5h5.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))
    # build a symmetric U column (sum of bond primitives only)
    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "bond":
            U[i, 0] = 1.0
    U_sym, info = symmetrize_u(
        U,
        prims,
        Z,
        coords,
        tol=1.0e-2,
        max_n=8,
        keep_a1_only=False,
        label_symmetry=True,
    )
    labels = info.get("labels", [])
    assert any(lab.get("label") == "A1" for lab in labels)


def test_a1_labels_ring_angles():
    # symmetric sum of ring angles should be A1-like
    xyz = Path(__file__).resolve().parent / "data" / "c5h5.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))

    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "angle":
            U[i, 0] = 1.0

    U_sym, info = symmetrize_u(
        U,
        prims,
        Z,
        coords,
        tol=1.0e-2,
        max_n=8,
        keep_a1_only=False,
        label_symmetry=True,
    )
    labels = info.get("labels", [])
    assert any(lab.get("label") == "A1" for lab in labels)


def test_a1_labels_naphthalene_angles_dihedrals():
    xyz = Path(__file__).resolve().parent / "data" / "naphthalene_c10.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))

    # A1 from sum of all valence angles
    U_ang = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "angle":
            U_ang[i, 0] = 1.0
    _, info_ang = symmetrize_u(
        U_ang,
        prims,
        Z,
        coords,
        tol=1.0e-3,
        max_n=8,
        keep_a1_only=False,
        label_symmetry=True,
    )
    labels_ang = info_ang.get("labels", [])
    assert any(lab.get("label") == "A1" for lab in labels_ang)

    # A1 from sum of all dihedrals (intra-cycle in planar system)
    U_dih = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "dihedral":
            U_dih[i, 0] = 1.0
    _, info_dih = symmetrize_u(
        U_dih,
        prims,
        Z,
        coords,
        tol=1.0e-3,
        max_n=8,
        keep_a1_only=False,
        label_symmetry=True,
    )
    labels_dih = info_dih.get("labels", [])
    assert any(lab.get("label") == "A1" for lab in labels_dih)


def test_symmetry_tol_h_allows_perturbed_h():
    # CH4 with one H slightly perturbed: should recover T-like symmetry with tol_H
    xyz = Path(__file__).resolve().parent / "data" / "ch4.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    coords_ang = coords_ang.copy()
    coords_ang[1] += np.array([0.01, 0.0, 0.0])
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    prims = primitives_from_topology(coords_ang / 0.52917721092, Z, np.deg2rad(170.0))

    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "bond":
            U[i, 0] = 1.0

    U_sym, info = symmetrize_u(
        U,
        prims,
        Z,
        coords_ang / 0.52917721092,
        tol=1.0e-3,
        max_n=10,
        keep_a1_only=False,
        label_symmetry=True,
        tol_H=0.10,
        quasi_tol=0.10,
        heavy_only_orient=True,
    )
    assert info.get("quasi_point_group", info.get("point_group")) in ("Td", "T", "Th", "Cs")


def test_symmetry_center_idx_spherical_guess():
    # SF6 with explicit center index should yield Oh
    from pathlib import Path
    from survibfit.modify_geom import read_xyz
    from topology.elements import atomic_number

    xyz = Path(__file__).resolve().parent / "data" / "sf6.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))

    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "bond":
            U[i, 0] = 1.0

    _, info = symmetrize_u(
        U,
        prims,
        Z,
        coords,
        tol=1.0e-3,
        max_n=10,
        keep_a1_only=False,
        label_symmetry=True,
        center_idx=0,
    )
    assert info.get("point_group") in ("Oh", "O")


def test_max_dev_strict_filters_ops():
    xyz = Path(__file__).resolve().parent / "data" / "ch4.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))
    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "bond":
            U[i, 0] = 1.0
    _, info = symmetrize_u(
        U,
        prims,
        Z,
        coords,
        tol=1.0e-3,
        max_n=10,
        max_dev_strict=1e-6,
    )
    # strict filter should keep at least identity
    assert "E" in info.get("elements", [])


def test_symmetry_profile_and_confidence():
    from pathlib import Path
    from survibfit.modify_geom import read_xyz
    from topology.elements import atomic_number
    from survibfit.symmetry_global import symmetrize_u
    from survibfit.pipeline import build_topology

    xyz = Path(__file__).resolve().parent / "data" / "ch4.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))

    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "bond":
            U[i, 0] = 1.0

    _, info = symmetrize_u(
        U,
        prims,
        Z,
        coords,
        tol=1.0e-3,
        max_n=10,
        profile=True,
    )
    perf = info.get("perf", {})
    assert perf.get("ops_total", 0) >= perf.get("ops_kept", 0)

    from survibfit.transforms import build_u
    _, _, ringset = build_topology(coords, Z)
    U_build = build_u(
        prims,
        coords,
        Z=Z,
        ringset=ringset,
        symmetrize_global=True,
        symmetry_confidence=True,
    )
    assert U_build.shape[0] == len(prims)


def test_ignore_isotopes_effect():
    # CH4 with one H replaced by D (simulate isotope label)
    from pathlib import Path
    from survibfit.modify_geom import read_xyz
    from topology.elements import atomic_number
    from survibfit.symmetry_global import symmetrize_u

    xyz = Path(__file__).resolve().parent / "data" / "ch4.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    atoms = list(atoms)
    atoms[1] = "D"
    Z = np.array([atomic_number("H") if a == "D" else atomic_number(a) for a in atoms], dtype=int)
    coords = coords_ang / 0.52917721092
    prims = primitives_from_topology(coords, Z, np.deg2rad(170.0))
    U = np.zeros((len(prims), 1), dtype=float)
    for i, p in enumerate(prims):
        if p.kind == "bond":
            U[i, 0] = 1.0

    symbols = atoms
    _, info_iso = symmetrize_u(U, prims, Z, coords, ignore_isotopes=False, symbols_override=symbols)
    _, info_ign = symmetrize_u(U, prims, Z, coords, ignore_isotopes=True, symbols_override=symbols)
    assert info_ign.get("point_group") in ("Td", "T", "Th", "O", "Oh")
    # without isotope ignore, symmetry should drop
    assert info_iso.get("point_group") in ("C1", "Cs", "C2", "C2v", "C2h", "C3v", "C3h")


def test_point_group_ch4_t():
    xyz = Path(__file__).resolve().parent / "data" / "ch4.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    symbols = [a for a in atoms]
    coords_oriented = orient_coords(coords_ang, weights=Z)
    elements, classes, permutations = symmetry_elements_from_geometry(
        symbols, coords_oriented, tol=1.0e-3, max_n=10
    )
    from survibfit.symmetry_global import _group_label
    pg = _group_label(elements)
    assert pg in ("Td", "T")


def test_point_group_sf6_o():
    xyz = Path(__file__).resolve().parent / "data" / "sf6.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    symbols = [a for a in atoms]
    coords_oriented = orient_coords(coords_ang, weights=Z)
    elements, classes, permutations = symmetry_elements_from_geometry(
        symbols, coords_oriented, tol=1.0e-3, max_n=10
    )
    from survibfit.symmetry_global import _group_label
    pg = _group_label(elements)
    assert pg in ("Oh", "O")
