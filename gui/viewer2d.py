"""
viewer2d.py

2D molecular drawing for Merlino GUI.

DESIGN RULE (FINAL):
- viewer2d is a TOPOLOGICAL viewer, not a geometric one
- viewer2d works ONLY if xyzin contains a #SMILES section
- viewer2d reads ONLY from xyzin
- no in-memory structures are ever used
"""

from __future__ import annotations

# ==============================================================
# Optional RDKit dependency
# ==============================================================

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdDepictor
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from PIL import Image
from pathlib import Path


# ==============================================================
# Utilities
# ==============================================================

def _extract_smiles_from_xyzin(xyzin_path: Path) -> str | None:
    """
    Extract SMILES from the #SMILES section of xyzin.
    Returns None if the section does not exist.
    """
    try:
        with open(xyzin_path, "r") as f:
            lines = f.readlines()
    except Exception:
        return None

    in_smiles = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.upper().startswith("#"):
            in_smiles = (line.upper() == "#SMILES")
            continue
        if in_smiles:
            return line

    return None


# ==============================================================
# Core functionality
# ==============================================================

def image_from_xyzin(xyzin_path: str | Path, size=(160, 160)) -> Image.Image | None:
    """
    Generate a 2D depiction from xyzin.

    Rules:
    - xyzin MUST exist
    - #SMILES section MUST exist
    - if any condition is not met, returns None
    """
    if not RDKIT_AVAILABLE:
        return None

    xyzin_path = Path(xyzin_path)

    smiles = _extract_smiles_from_xyzin(xyzin_path)
    if not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None

        rdDepictor.Compute2DCoords(mol)

        img = Draw.MolToImage(
            mol,
            size=size,
            kekulize=True,
            wedgeBonds=True
        )
        return img if isinstance(img, Image.Image) else None

    except Exception:
        return None


# ==============================================================
# Validation helper (optional)
# ==============================================================

def has_smiles(xyzin_path: str | Path) -> bool:
    """
    Return True if xyzin contains a #SMILES section.
    """
    xyzin_path = Path(xyzin_path)
    return _extract_smiles_from_xyzin(xyzin_path) is not None
