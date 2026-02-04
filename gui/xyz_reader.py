"""
XYZ reader for Merlino 3.0.

Behavior:
- Reads an XYZ (standard or enriched)
- Replaces ONLY the XYZ block in working/xyzin
- Preserves all sections (#BASIC, #SMILES, ...)
- Does not return anything
"""

from pathlib import Path

from .xyzin_utils import replace_xyz_block


# --------------------------------------------------
# Atomic number -> symbol map (1..118)
# --------------------------------------------------
Z2SYM = {
    1: "H",   2: "He",
    3: "Li",  4: "Be",  5: "B",   6: "C",   7: "N",   8: "O",   9: "F",   10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",  17: "Cl", 18: "Ar",
    19: "K",  20: "Ca", 21: "Sc", 22: "Ti", 23: "V",  24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y",  40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru",
    45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd",
    49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I",  54: "Xe",
    55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm",
    63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu",
    72: "Hf", 73: "Ta", 74: "W",  75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au",
    80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn",
    87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",  93: "Np", 94: "Pu",
    95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm", 101: "Md",
    102: "No", 103: "Lr",
    104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt",
    110: "Ds", 111: "Rg", 112: "Cn",
    113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
}

VALID_SYMBOLS = set(Z2SYM.values())


def _normalize_symbol(sym: str) -> str:
    """
    Normalize element symbol capitalization:
      'c'  -> 'C'
      'cl' -> 'Cl'
      'CL' -> 'Cl'
      'cL' -> 'Cl'
    """
    s = sym.strip()
    if not s:
        raise ValueError("Empty atomic symbol in XYZ")

    s = s.lower()
    if len(s) == 1:
        return s.upper()
    return s[0].upper() + s[1:]


def _convert_atom_token_to_symbol(token: str) -> str:
    """
    Convert first token of an XYZ coordinate line to a valid atomic symbol.

    Accepts:
    - Symbol (any capitalization): 'C', 'cl', 'SI', ...
    - Atomic number: '6', '17', ...

    Returns:
    - Normalized symbol: 'C', 'Cl', 'Si', ...

    Raises:
    - ValueError if the symbol/atomic number is not recognized.
    """
    tok = token.strip()
    if not tok:
        raise ValueError("Missing element field in XYZ line")

    # Atomic number
    if tok.isdigit():
        z = int(tok)
        if z not in Z2SYM:
            raise ValueError(f"Invalid atomic number in XYZ: {z}")
        return Z2SYM[z]

    # Symbol (case-insensitive)
    sym = _normalize_symbol(tok)
    if sym not in VALID_SYMBOLS:
        raise ValueError(f"Invalid atomic symbol in XYZ: '{tok}' (normalized '{sym}')")
    return sym


def read_xyz(path):
    """
    Read an XYZ file and replace the XYZ block in xyzin.

    Parameters
    ----------
    path : str or Path
        Path to XYZ file
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    # --------------------------------------------------
    # IMPORTANT: DO NOT remove empty lines!
    # XYZ standard allows an empty comment line (line 2).
    # --------------------------------------------------
    raw_lines = path.read_text().splitlines()

    if len(raw_lines) < 2:
        raise ValueError("Invalid XYZ file: too few lines")

    # --------------------------------------------------
    # Atom count (1st line)
    # --------------------------------------------------
    try:
        nat = int(raw_lines[0].strip())
    except ValueError:
        raise ValueError("Invalid XYZ file: first line is not atom count")

    # --------------------------------------------------
    # Comment line (2nd line) - can be empty
    # --------------------------------------------------
    comment = raw_lines[1].rstrip()

    # --------------------------------------------------
    # Parse coordinate lines starting from line 3
    # Skip blank lines, but DO NOT treat them as data
    # --------------------------------------------------
    coord_lines = []
    for line in raw_lines[2:]:
        if not line.strip():
            continue

        tokens = line.split()
        if len(tokens) < 4:
            continue

        sym = _convert_atom_token_to_symbol(tokens[0])

        # rebuild line with normalized symbol + original rest
        new_line = " ".join([sym] + tokens[1:])
        coord_lines.append(new_line)

        if len(coord_lines) == nat:
            break

    if len(coord_lines) != nat:
        raise ValueError(
            f"Invalid XYZ file: expected {nat} coordinate lines, found {len(coord_lines)}"
        )

    xyz_block = [str(nat), comment] + coord_lines

    # --------------------------------------------------
    # Update xyzin
    # --------------------------------------------------
    replace_xyz_block(xyz_block)

