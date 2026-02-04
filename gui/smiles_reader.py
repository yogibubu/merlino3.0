"""
SMILES reader for Merlino 3.0.

Behavior:
- Uses RDKit if available
- Updates XYZ block
- Updates #SMILES section
- Does NOT touch #BASIC
- Raises exceptions on invalid input
"""

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False

from .xyzin_utils import replace_xyz_block, append_section


def read_smiles(smiles: str, **_):
    """
    Read a SMILES string and update xyzin.

    Parameters
    ----------
    smiles : str
        SMILES string (must be valid)
    """
    if not smiles or not isinstance(smiles, str):
        raise ValueError("Empty SMILES")

    smiles = smiles.strip()

    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit not available")

    # --------------------------------------------------
    # Parse SMILES
    # --------------------------------------------------
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            raise ValueError("Cannot parse SMILES")

        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        mol = Chem.AddHs(mol)
    except Exception as err:
        raise ValueError(f"Invalid SMILES: {err}")

    # --------------------------------------------------
    # 3D embedding
    # --------------------------------------------------
    try:
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        pass

    # --------------------------------------------------
    # Build XYZ block
    # --------------------------------------------------
    try:
        conf = mol.GetConformer()
        atoms = mol.GetAtoms()

        xyz = [str(len(atoms)), "Generated from SMILES"]
        for i, atom in enumerate(atoms):
            pos = conf.GetAtomPosition(i)
            xyz.append(
                f"{atom.GetSymbol():2s} "
                f"{pos.x:15.8f} {pos.y:15.8f} {pos.z:15.8f}"
            )
    except Exception as err:
        raise RuntimeError(f"Cannot build XYZ: {err}")

    # --------------------------------------------------
    # Write to xyzin (canonical behavior)
    # --------------------------------------------------
    replace_xyz_block(xyz)
    append_section("SMILES", [smiles])

