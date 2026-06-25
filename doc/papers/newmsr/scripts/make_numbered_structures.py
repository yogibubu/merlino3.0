from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/Users/vincenzobarone/merlino3.0")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from merlino_fit.topology.elements import atomic_number
from merlino_fit.topology.pipeline import build_topology_objects


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "figures" / "data"
OUT = ROOT / "figures"


COLORS = {
    "C": "#2f3437",
    "H": "#f7f3ea",
    "N": "#2b67b1",
    "O": "#c74432",
}


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    nat = int(lines[0].strip())
    atoms: list[str] = []
    coords: list[list[float]] = []
    for raw in lines[2 : 2 + nat]:
        parts = raw.split()
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return atoms, np.asarray(coords, dtype=float)


def topology_bonds(atoms: list[str], coords: np.ndarray) -> list[tuple[int, int]]:
    z_numbers = np.array([atomic_number(atom) for atom in atoms], dtype=int)
    _continuous, graph, _ringset, _synthons, _aromaticity = build_topology_objects(coords, z_numbers)
    return sorted(tuple(sorted(pair)) for pair in graph.bonds)


def project(coords: np.ndarray) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    _u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    xy = centered @ vt[:2].T
    if len(singular) > 2 and singular[2] > 1.0e-3:
        # Separate out-of-plane atoms in near-planar molecules without changing numbering.
        xy[:, 0] += 1.10 * (centered @ vt[2])
    if np.ptp(xy[:, 0]) < np.ptp(xy[:, 1]):
        xy = xy[:, ::-1]
    return xy


def draw_structure(name: str, xyz: Path, output: Path) -> None:
    atoms, coords = read_xyz(xyz)
    xy = project(coords)
    bonds = topology_bonds(atoms, coords)
    manual_xy = manual_layout(name, atoms)
    if manual_xy is not None:
        xy = manual_xy
    display_xy = xy.copy()

    fig, ax = plt.subplots(figsize=(4.8, 4.1))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor("white")

    span = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1])), 1.0)
    min_h_distance = 0.18 * span
    adjacency = {idx: [] for idx in range(len(atoms))}
    for i, j in bonds:
        adjacency[i].append(j)
        adjacency[j].append(i)
    for idx, atom in enumerate(atoms):
        if atom != "H" or len(adjacency[idx]) != 1:
            continue
        parent = adjacency[idx][0]
        direction = display_xy[idx] - display_xy[parent]
        norm = float(np.linalg.norm(direction))
        if norm >= min_h_distance:
            continue
        if norm < 1.0e-8:
            direction = display_xy[parent] - display_xy.mean(axis=0)
            norm = float(np.linalg.norm(direction))
        if norm < 1.0e-8:
            direction = np.array([1.0, 0.0])
            norm = 1.0
        display_xy[idx] = display_xy[parent] + direction / norm * min_h_distance

    for i, j in bonds:
        ax.plot(
            [display_xy[i, 0], display_xy[j, 0]],
            [display_xy[i, 1], display_xy[j, 1]],
            color="#262626",
            lw=1.8,
            solid_capstyle="round",
            zorder=1,
        )

    atom_radius = 0.070 * span
    label_radius = 0.095 * span

    for idx, (atom, pos) in enumerate(zip(atoms, display_xy), start=1):
        face = COLORS.get(atom, "#d9d9d9")
        edge = "#262626"
        radius = atom_radius * (0.78 if atom == "H" else 1.0)
        circle = plt.Circle(pos, radius, facecolor=face, edgecolor=edge, lw=1.1, zorder=3)
        ax.add_patch(circle)
        text_color = "white" if atom != "H" else "#222222"
        ax.text(
            pos[0],
            pos[1],
            f"{atom}{idx}",
            ha="center",
            va="center",
            fontsize=8.2 if atom == "H" else 8.8,
            color=text_color,
            fontweight="bold",
            zorder=4,
        )

    margin = label_radius * 1.6
    ax.set_xlim(float(display_xy[:, 0].min() - margin), float(display_xy[:, 0].max() + margin))
    ax.set_ylim(float(display_xy[:, 1].min() - margin), float(display_xy[:, 1].max() + margin))
    ax.set_title(name, fontsize=11, pad=8)
    fig.savefig(output, bbox_inches="tight", transparent=True)
    plt.close(fig)


def manual_layout(name: str, atoms: list[str]) -> np.ndarray | None:
    if name == "Norcamphor" and len(atoms) == 18:
        return np.asarray(
            [
                [0.00, 1.35],   # C1
                [-1.05, 0.70],  # C2
                [-0.85, -0.25], # C3
                [0.20, -0.80],  # C4
                [1.20, -0.25],  # C5
                [1.05, 0.90],   # C6
                [0.25, 0.20],   # C7
                [-1.95, 1.00],  # O8
                [-0.10, 2.15],  # H9
                [-1.60, -0.55], # H10
                [-0.90, -1.10], # H11
                [0.20, -1.55],  # H12
                [-0.45, 0.10],  # H13
                [0.55, 0.75],   # H14
                [1.85, -0.55],  # H15
                [1.25, -1.05],  # H16
                [1.65, 1.35],   # H17
                [1.85, 0.75],   # H18
            ],
            dtype=float,
        )
    if name == "Succinic anhydride" and len(atoms) == 11:
        return np.asarray(
            [
                [0.00, -1.25],   # O1
                [1.15, -0.45],   # C2
                [0.70, 0.92],    # C3
                [-0.70, 0.92],   # C4
                [-1.15, -0.45],  # C5
                [2.15, -0.05],   # O6
                [-2.15, -0.05],  # O7
                [1.15, 1.78],    # H8
                [0.34, 2.12],    # H9
                [-1.15, 1.78],   # H10
                [-0.34, 2.12],   # H11
            ],
            dtype=float,
        )
    return None


def main() -> None:
    draw_structure(
        "Glycolaldehyde",
        DATA / "glycolaldehyde_parent.xyz",
        OUT / "glycolaldehyde_numbering.pdf",
    )
    draw_structure(
        "Glycine I",
        DATA / "glycine_I_parent.xyz",
        OUT / "glycine_I_numbering.pdf",
    )
    draw_structure(
        "Glycine II",
        DATA / "glycine_II_parent.xyz",
        OUT / "glycine_II_numbering.pdf",
    )
    draw_structure(
        "Cyclopentadiene",
        DATA / "cyclopentadiene_parent.xyz",
        OUT / "cyclopentadiene_numbering.pdf",
    )
    draw_structure(
        "Nitrobenzene",
        DATA / "nitrobenzene_parent.xyz",
        OUT / "nitrobenzene_numbering.pdf",
    )
    draw_structure(
        "p-EBN",
        DATA / "p-EBN_parent.xyz",
        OUT / "p-EBN_numbering.pdf",
    )
    draw_structure(
        "Azulene",
        DATA / "azulene_parent.xyz",
        OUT / "azulene_numbering.pdf",
    )
    draw_structure(
        "Norcamphor",
        DATA / "norcamphor_parent.xyz",
        OUT / "norcamphor_numbering.pdf",
    )
    draw_structure(
        "Maleic anhydride",
        DATA / "maleic_anhydride_parent.xyz",
        OUT / "maleic_anhydride_numbering.pdf",
    )
    draw_structure(
        "Phthalic anhydride",
        DATA / "phthalic_anhydride_parent.xyz",
        OUT / "phthalic_anhydride_numbering.pdf",
    )
    draw_structure(
        "Succinic anhydride",
        DATA / "succinic_anhydride_parent.xyz",
        OUT / "succinic_anhydride_numbering.pdf",
    )


if __name__ == "__main__":
    main()
