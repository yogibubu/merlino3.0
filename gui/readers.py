"""
Reader dispatcher for Merlino 3.0.

FINAL CONTRACT
--------------
- This module dispatches input sources to the correct reader
- Readers:
    - write ONLY to working/xyzin
    - return None
- No path logic
- No chemistry logic
"""

from .smiles_reader import read_smiles
from .xyz_reader import read_xyz
from .gaussian import read_gaussian
from .molpro import read_molpro
from .mrcc import read_mrcc

def read_structure(input_source):
    """
    Dispatch input source to the appropriate reader.

    Parameters
    ----------
    input_source : object
        Must have attribute 'kind'
        Additional attributes depend on kind:
            - smiles : str
            - path   : str or Path
            - strict : bool (optional, for SMILES)
    """

    if not hasattr(input_source, "kind"):
        raise TypeError("input_source must have attribute 'kind'")

    kind = input_source.kind.lower()

    # --------------------------------------------------
    # SMILES
    # --------------------------------------------------
    if kind == "smiles":
        if not hasattr(input_source, "smiles"):
            raise AttributeError("SMILES input requires 'smiles' attribute")

        strict = getattr(input_source, "strict", False)
        read_smiles(input_source.smiles, strict=strict)
        return

    # --------------------------------------------------
    # XYZ
    # --------------------------------------------------
    if kind == "xyz":
        if not hasattr(input_source, "path"):
            raise AttributeError("XYZ input requires 'path' attribute")
        read_xyz(input_source.path)
        return

    # --------------------------------------------------
    # GAUSSIAN
    # --------------------------------------------------
    if kind == "gaussian":
        if not hasattr(input_source, "path"):
            raise AttributeError("Gaussian input requires 'path' attribute")
        read_gaussian(input_source.path)
        return

    # --------------------------------------------------
    # MOLPRO
    # --------------------------------------------------
    if kind == "molpro":
        if not hasattr(input_source, "path"):
            raise AttributeError("Molpro input requires 'path' attribute")
        read_molpro(input_source.path)
        return

    # --------------------------------------------------
    # MRCC
    # --------------------------------------------------
    if kind == "mrcc":
        if not hasattr(input_source, "path"):
            raise AttributeError("MRCC input requires 'path' attribute")
        read_mrcc(input_source.path)
        return

    # --------------------------------------------------
    # Unsupported
    # --------------------------------------------------
    raise ValueError(f"Unsupported input kind: {input_source.kind}")
