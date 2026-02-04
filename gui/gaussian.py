"""
Gaussian output reader (xyzin-only).

FINAL CONTRACT (Merlino 3.0):
- This reader WRITES xyzin
- NO Structure objects
- NO return values
 - Parsing preserved

PATCH:
- If input is .fchk/.fch, COPY it into working/fchkin so downstream modules
  (rotational/vibrational/qcent/...) can use it.
- If input is .log/.out, COPY it into working/logfin and always parse from
  the copied files.
"""

from pathlib import Path
from typing import List, Tuple
import shutil
import re

from .xyzin_utils import replace_xyz_block, append_section, read_xyzin_lines, write_xyzin_lines, remove_section

BOHR_TO_ANG = 0.52917721092
# Keep consistent with topology/discrete_graph.py
BOND_ORDER_WRITE_THRESHOLD = 0.2


# ==================================================
# Periodic table (Z -> symbol)
# ==================================================
def _atomic_symbol(z: int) -> str:
    periodic = [
        "",
        "H", "He",
        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
        "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn",
    ]

    if z < 1 or z >= len(periodic):
        raise ValueError(f"Unsupported atomic number: {z}")

    return periodic[z]


# ==================================================
# Parsing helpers (UNCHANGED)
# ==================================================
def _parse_charge_mult(lines: List[str]) -> Tuple[int, int]:
    for line in lines:
        if "Charge =" in line and "Multiplicity =" in line:
            parts = line.replace("=", "").split()
            charge = int(parts[parts.index("Charge") + 1])
            mult = int(parts[parts.index("Multiplicity") + 1])
            return charge, mult
    return 0, 1


def _parse_orientations(lines: List[str]):
    standard_geoms = []
    input_geoms = []

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        is_standard = "Standard orientation:" in line
        is_input = "Input orientation:" in line

        if is_standard or is_input:
            i += 1
            while i < n and not lines[i].strip().startswith("----"):
                i += 1
            i += 1
            while i < n and not lines[i].strip().startswith("----"):
                i += 1
            i += 1

            geom = []

            while i < n:
                l = lines[i].strip()
                if not l or l.startswith("----"):
                    break

                fields = l.split()
                if len(fields) < 6:
                    break

                atomic_number = int(fields[1])
                x, y, z = map(float, fields[3:6])
                symbol = _atomic_symbol(atomic_number)

                geom.append((symbol, x, y, z))
                i += 1

            if geom:
                if is_standard:
                    standard_geoms.append(geom)
                else:
                    input_geoms.append(geom)

        i += 1

    return standard_geoms if standard_geoms else input_geoms


# ==================================================
# Public API
# ==================================================

def _find_gaussian_companion_files(path: Path):
    """Try to locate Gaussian LOG/OUT and FCHK/FCH with the same basename in the same folder."""
    path = Path(path)
    ext = path.suffix.lower()
    stem = path.stem
    folder = path.parent

    log_path = None
    fchk_path = None

    if ext in {".log", ".out"}:
        log_path = path
        cand1 = folder / (stem + ".fchk")
        cand2 = folder / (stem + ".fch")
        if cand1.exists():
            fchk_path = cand1
        elif cand2.exists():
            fchk_path = cand2

    elif ext in {".fchk", ".fch"}:
        fchk_path = path
        cand1 = folder / (stem + ".log")
        cand2 = folder / (stem + ".out")
        if cand1.exists():
            log_path = cand1
        elif cand2.exists():
            log_path = cand2

    return log_path, fchk_path


def _parse_deltabvib_from_log(log_path: Path):
    """Parse the exact anharmonic rotational constants block (in MHz) and compute DeltaBvib."""
    try:
        lines = Path(log_path).read_text(errors="ignore").splitlines()
    except Exception:
        return None, "log read error"

    # Find the exact header line
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Rotational Constants (in MHz)":
            idx = i
            break
    if idx is None:
        return None, "Gaussian calculation is not anharmonic (Rotational Constants in MHz block not found)"

    # We expect:
    # idx+1: dashes
    # idx+2: description line
    # idx+3..5: Ae/Be/Ce lines (order fixed)
    def _parse_pair(line, left_key, right_key):
        # Example: 'Ae=  698577.947     A00=  707186.293     A0=  707059.182'
        try:
            # Split by whitespace, keep tokens like 'Ae=' and numbers.
            toks = line.replace("=", "= ").split()
            # Build dict key->value for Ae and A00
            d = {}
            for j in range(len(toks) - 1):
                if toks[j].endswith("=") and toks[j+1].replace(".","",1).replace("-","",1).isdigit() is False:
                    pass
            # Simpler deterministic scan:
            val_left = None
            val_right = None
            for j, t in enumerate(toks):
                if t == left_key + "=" and j + 1 < len(toks):
                    val_left = float(toks[j + 1])
                if t == right_key + "=" and j + 1 < len(toks):
                    val_right = float(toks[j + 1])
            if val_left is None or val_right is None:
                return None, None
            return val_left, val_right
        except Exception:
            return None, None

    # Search forward for Ae line after header within next ~20 lines
    Ae = Be = Ce = A00 = B00 = C00 = None
    for k in range(idx, min(idx + 30, len(lines))):
        line = lines[k].strip()
        if line.startswith("Ae="):
            Ae, A00 = _parse_pair(line, "Ae", "A00")
        elif line.startswith("Be="):
            Be, B00 = _parse_pair(line, "Be", "B00")
        elif line.startswith("Ce="):
            Ce, C00 = _parse_pair(line, "Ce", "C00")
        if Ae is not None and Be is not None and Ce is not None and A00 is not None and B00 is not None and C00 is not None:
            break

    if Ae is None or Be is None or Ce is None or A00 is None or B00 is None or C00 is None:
        return None, "Rotational Constants (in MHz) block found, but parsing failed"

    # Values in the log block are in MHz.
    dva = (A00 - Ae) 
    dvb = (B00 - Be) 
    dvc = (C00 - Ce) 
    return (dva, dvb, dvc), None


def _has_output_pickett(lines: List[str]) -> bool:
    for line in lines:
        if "output=pickett" in line.lower():
            return True
    return False


def _parse_gaussian_float(token: str):
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except Exception:
        return None


def _parse_harmonic_freq_ir_from_log(lines: List[str]):
    """
    Parse harmonic frequencies (cm^-1) and IR intensities (KM/Mole) from Gaussian log.
    Uses the LAST "Harmonic frequencies (cm**-1)" block.
    Returns (freq_cm1, ir_km_mol) lists or (None, None) if not found.
    """
    header_idx = None
    for i, line in enumerate(lines):
        if "Harmonic frequencies (cm**-1)" in line:
            header_idx = i

    if header_idx is None:
        return None, None

    freq = []
    ir = []

    for k in range(header_idx + 1, len(lines)):
        s = lines[k].strip()
        if not s:
            continue

        if s.startswith("----") or s.startswith("Thermochemistry") or s.startswith("Temperature"):
            break

        if "Frequencies --" in s:
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[DEde][+-]?\d+)?", s)
            for n in nums:
                v = _parse_gaussian_float(n)
                if v is not None:
                    freq.append(v)
            continue

        if s.startswith("IR Inten"):
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[DEde][+-]?\d+)?", s)
            for n in nums:
                v = _parse_gaussian_float(n)
                if v is not None:
                    ir.append(v)
            continue

        # Stop once we enter the normal coordinates table
        if s.startswith("Atom  AN"):
            continue

    if not freq:
        return None, None

    return freq, (ir if ir else None)


def _parse_anharmonic_x_matrix(lines: List[str]):
    """
    Parse the "Total Anharmonic X Matrix (in cm^-1)" block.
    Returns full symmetric matrix as list of lists, or None if not found.
    """
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Total Anharmonic X Matrix (in cm^-1)":
            header_idx = i

    if header_idx is None:
        return None

    k = header_idx + 1
    while k < len(lines) and (not lines[k].strip() or lines[k].strip().startswith("-")):
        k += 1
    if k >= len(lines):
        return None

    col_line = lines[k].strip()
    cols = [int(x) for x in re.findall(r"\d+", col_line)]
    if not cols:
        return None
    n = max(cols)

    chi = [[None for _ in range(n)] for _ in range(n)]

    for i in range(k + 1, len(lines)):
        s = lines[i].strip()
        if not s or s.startswith("=") or s.startswith("-"):
            break
        parts = s.split()
        if not parts or not parts[0].isdigit():
            break
        row = int(parts[0])
        vals = parts[1:]
        for j, tok in enumerate(vals, start=1):
            val = _parse_gaussian_float(tok)
            if val is None:
                continue
            if j <= n and row <= n:
                chi[row - 1][j - 1] = val
                chi[j - 1][row - 1] = val

    # If we parsed nothing, return None
    any_val = any(chi[i][j] is not None for i in range(n) for j in range(n))
    if not any_val:
        return None

    # Fill missing diagonal with 0.0 if absent
    for i in range(n):
        if chi[i][i] is None:
            chi[i][i] = 0.0

    return chi


def _build_vibrational_lines_from_log(lines: List[str]):
    """
    Build #VIBRATIONAL section lines from Gaussian LOG.
    Includes:
      - harmonic frequencies (cm^-1)
      - harmonic IR intensities (KM/Mole)
      - anharmonic X (chi) matrix (cm^-1), if present
    """
    freq, ir = _parse_harmonic_freq_ir_from_log(lines)
    if not freq:
        return None

    out = []
    out.append("freq_cm1 = " + " ".join(f"{x:.6f}" for x in freq))

    if ir is not None and len(ir) == len(freq):
        out.append("ir_inten_km_mol = " + " ".join(f"{x:.6f}" for x in ir))

    chi = _parse_anharmonic_x_matrix(lines)
    if chi is not None:
        out.append("chi_cm1 = [")
        n = len(chi)
        for i in range(n):
            for j in range(i + 1):
                if chi[i][j] is None:
                    continue
                out.append(f"{i+1:3d} {j+1:3d} {chi[i][j]: .8f}")
        out.append("]")

    return out if out else None


def _parse_rot_constants_mhz(lines: List[str]):
    """
    Parse "Rotational Constants (in MHz)" block.
    Returns dict like {"Ae": val, "A00": val, "A0": val, ...}
    """
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Rotational Constants (in MHz)":
            idx = i
            break
    if idx is None:
        return {}

    out = {}
    pattern = re.compile(r"([A-Z][A-Za-z0-9]*)=\s*([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)")
    for k in range(idx, min(idx + 15, len(lines))):
        line = lines[k].strip()
        if not line or line.startswith("----"):
            continue
        if line.startswith(("Ae=", "Be=", "Ce=")):
            for m in pattern.finditer(line):
                key = m.group(1)
                val = _parse_gaussian_float(m.group(2))
                if val is not None:
                    out[key] = val
    return out


def _parse_spfit_like_block(lines: List[str]):
    """
    Parse a SPFIT-like block if present.
    Returns dict with keys like A,B,C,DJ,DJK,DK,d1,d2,DelJ,DelJK,DelK,delJ,delK
    """
    idx = None
    for i, line in enumerate(lines):
        if "SPFIT-like block" in line:
            idx = i
            break
    if idx is None:
        return {}

    out = {}
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9]*)\s*=\s*([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)")
    for k in range(idx + 1, min(idx + 60, len(lines))):
        line = lines[k].strip()
        if not line:
            continue
        m = pattern.match(line)
        if not m:
            continue
        key = m.group(1)
        val = _parse_gaussian_float(m.group(2))
        if val is not None:
            out[key] = val
    return out


def _parse_reduced_constants(lines: List[str], header: str):
    """
    Parse Gaussian reduced Hamiltonian constants (MHz).
    header: exact header line to locate.
    Returns dict with canonical keys (DelJ_MHz, DJ_MHz, etc).
    """
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == header:
            idx = i
            break
    if idx is None:
        return {}

    out = {}
    number_pattern = re.compile(r"[+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?")
    for k in range(idx + 1, min(idx + 20, len(lines))):
        raw = lines[k].rstrip()
        if not raw.strip():
            break
        if ":" not in raw:
            continue
        label = raw.split(":", 1)[0].strip()
        nums = number_pattern.findall(raw)
        if not nums:
            continue
        val = _parse_gaussian_float(nums[-1])
        if val is None:
            continue

        label_clean = re.sub(r"\s+", "", label)

        if label_clean.startswith("DELTA"):
            key = {
                "DELTAJ": "DelJ_MHz",
                "DELTAK": "DelK_MHz",
                "DELTAJK": "DelJK_MHz",
            }.get(label_clean.upper())
        elif label_clean.lower().startswith("delta"):
            key = {
                "deltaj": "delJ_MHz",
                "deltak": "delK_MHz",
                "deltajk": "delJK_MHz",
            }.get(label_clean.lower())
        else:
            key = {
                "dj": "DJ_MHz",
                "djk": "DJK_MHz",
                "dk": "DK_MHz",
                "d1": "d1_MHz",
                "d2": "d2_MHz",
            }.get(label_clean.lower())

        if key:
            out[key] = val

    return out


def _parse_pickett_dipole_debye(lines: List[str]):
    """
    Parse the Pickett-style "Dipole moment (Debye)" block.
    Returns (a, b, c) components if available.
    """
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Dipole moment (Debye):":
            idx = i
    if idx is None:
        return None

    number_pattern = re.compile(r"[+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?")
    for k in range(idx + 1, min(idx + 6, len(lines))):
        raw = lines[k].strip()
        if not raw:
            continue
        nums = number_pattern.findall(raw)
        if len(nums) >= 3:
            a = _parse_gaussian_float(nums[0])
            b = _parse_gaussian_float(nums[1])
            c = _parse_gaussian_float(nums[2])
            if a is not None and b is not None and c is not None:
                return (a, b, c)
    return None


def _parse_pickett_quadrupole_chi(lines: List[str]):
    """
    Parse Pickett-style nuclear quadrupole coupling constants [Chi] (MHz).
    Returns dict of Chi_* keys if available.
    """
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Nuclear quadrupole coupling constants [Chi] (MHz):":
            idx = i
    if idx is None:
        return {}

    out = {}
    number_pattern = re.compile(r"[+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?")
    kv_pattern = re.compile(r"(Chi\w+)\s*=?\s*([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)", re.I)

    for k in range(idx + 1, min(idx + 25, len(lines))):
        raw = lines[k].strip()
        if not raw:
            continue
        if raw.startswith(("Dipole moment", "Quartic", "Sextic", "Atoms with significant")):
            break

        matches = kv_pattern.findall(raw)
        if matches:
            for key, val in matches:
                v = _parse_gaussian_float(val)
                if v is not None:
                    out[f"{key}_MHz"] = v
            continue

        nums = number_pattern.findall(raw)
        if len(nums) >= 3:
            toks = raw.split()
            idx_tok = toks[0] if toks and toks[0].isdigit() else None
            elem_tok = toks[1] if len(toks) > 1 and toks[1].isalpha() else None
            if idx_tok:
                prefix = f"Chi{idx_tok}"
                if elem_tok:
                    prefix += f"_{elem_tok}"
            elif elem_tok:
                prefix = f"Chi_{elem_tok}"
            else:
                prefix = "Chi"

            a = _parse_gaussian_float(nums[0])
            b = _parse_gaussian_float(nums[1])
            c = _parse_gaussian_float(nums[2])
            if a is not None and b is not None and c is not None:
                out[f"{prefix}_aa_MHz"] = a
                out[f"{prefix}_bb_MHz"] = b
                out[f"{prefix}_cc_MHz"] = c

    return out


def _parse_pickett_sextic_constants(lines: List[str]):
    """
    Parse Pickett-style sextic centrifugal distortion constants (MHz).
    Returns dict of Phi*/H*/phi*/h* keys.
    """
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == "Sextic Centrifugal Distortion Constants (MHz)":
            idx = i
    if idx is None:
        return {}

    out = {}
    pattern = re.compile(r"([A-Za-z]+)\s*([A-Za-z0-9]+)?\s*=\s*([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)")

    for k in range(idx + 1, min(idx + 20, len(lines))):
        raw = lines[k].strip()
        if not raw:
            continue
        for m in pattern.finditer(raw):
            p1 = m.group(1)
            p2 = m.group(2) or ""
            label = f"{p1}{p2}".replace(" ", "")
            val = _parse_gaussian_float(m.group(3))
            if val is None:
                continue

            key = None
            if label.startswith("phi"):
                key = f"phi{label[3:]}_MHz"
            elif label.startswith("h"):
                key = f"h{label[1:]}_MHz"
            elif label.upper().startswith("PHI"):
                key = f"Phi{label[3:]}_MHz"
            elif label.upper().startswith("H"):
                key = f"H{label[1:]}_MHz"

            if key:
                out[key] = val

    return out


def _parse_cm5_charges(lines: List[str]):
    """
    Parse CM5 atomic charges from Gaussian log when present.
    Returns mapping {atom_index_1based: cm5_charge}.
    """
    header_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "cm5 charges" in low:
            header_idx = i

    if header_idx is None:
        return {}

    out = {}
    row_pat = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\b")
    num_pat = re.compile(r"[+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?")
    started = False

    for k in range(header_idx + 1, min(header_idx + 180, len(lines))):
        s = lines[k].strip()
        if not s:
            if started:
                break
            continue
        if s.startswith(("----", "Tot", "Sum")):
            if started:
                break
            continue

        m = row_pat.match(s)
        if not m:
            if started and s.lower().startswith(("hirshfeld", "mulliken", "natural")):
                break
            continue

        idx = int(m.group(1))
        nums = num_pat.findall(s[m.end():])
        if not nums:
            continue
        v = _parse_gaussian_float(nums[-1])
        if v is None:
            continue
        out[idx] = v
        started = True

    return out


def _parse_bond_order_pairs_from_line(line: str):
    """
    Parse bond-order pair entries from one line.
    Supports common Gaussian pair formats such as:
      B( 1-C, 2-C)=0.98
      B(1,2)=0.98
    """
    out = []
    patterns = [
        re.compile(
            r"B\(\s*(\d+)\s*-[A-Za-z]{1,2}\s*,\s*(\d+)\s*-[A-Za-z]{1,2}\s*\)\s*=\s*([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)"
        ),
        re.compile(
            r"B\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)"
        ),
    ]
    for pat in patterns:
        for m in pat.finditer(line):
            i = int(m.group(1))
            j = int(m.group(2))
            v = _parse_gaussian_float(m.group(3))
            if v is None:
                continue
            out.append((i, j, v))
    return out


def _parse_mayer_bond_orders(lines: List[str]):
    """
    Parse Mayer bond orders when present.
    Returns mapping {(i,j)_1based_sorted: bo}.
    """
    idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "mayer bond orders" in low or "mayer atomic bond orders" in low:
            idx = i

    if idx is None:
        return {}

    out = {}
    num_pat = re.compile(r"[+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?")
    cols_pat = re.compile(r"^\s*(\d+(?:\s+\d+)*)\s*$")
    row_pat = re.compile(r"^\s*(\d+)\s+[A-Za-z]{1,3}\s+(.+?)\s*$")
    current_cols = None
    parsed_any_matrix = False

    for k in range(idx + 1, min(idx + 260, len(lines))):
        s = lines[k]
        if not s.strip():
            continue
        low = s.lower()
        if low.startswith(("wiberg", "natural", "mulliken", "hirshfeld")) and out:
            break

        # Pair form (e.g., B( 1-C, 2-H)=0.9123)
        for i1, j1, v in _parse_bond_order_pairs_from_line(s):
            key = (i1, j1) if i1 < j1 else (j1, i1)
            out[key] = v

        # Matrix form (Atomic Valencies and Mayer Atomic Bond Orders)
        mcols = cols_pat.match(s)
        if mcols:
            current_cols = [int(x) for x in mcols.group(1).split()]
            parsed_any_matrix = True
            continue

        mrow = row_pat.match(s)
        if mrow and current_cols:
            i1 = int(mrow.group(1))
            nums = num_pat.findall(mrow.group(2))
            for j1, tok in zip(current_cols, nums):
                v = _parse_gaussian_float(tok)
                if v is None or i1 == j1:
                    continue
                key = (i1, j1) if i1 < j1 else (j1, i1)
                out[key] = v
            continue

        if parsed_any_matrix and out:
            # Stop when we leave the matrix block.
            if not low.startswith(("----",)):
                break
    return out


def _parse_total_bond_orders(lines: List[str]):
    """
    Parse total bond orders between atoms (if present in pair form).
    Returns mapping {(i,j)_1based_sorted: bo}.
    """
    idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "total bond order" in low and "atom" in low:
            idx = i

    if idx is None:
        return {}

    out = {}
    for k in range(idx + 1, min(idx + 260, len(lines))):
        s = lines[k]
        if not s.strip():
            continue
        if s.strip().startswith(("----", "Normal termination")) and out:
            break
        for i1, j1, v in _parse_bond_order_pairs_from_line(s):
            key = (i1, j1) if i1 < j1 else (j1, i1)
            out[key] = v
    return out


def _build_gaussian_topology_lines(lines: List[str]):
    """
    Build #GAUSSIAN_TOPOLOGY lines with optional CM5 charges and bond orders.
    Bond-order priority: Mayer > Total.
    Only bonded pairs are written (BO >= BOND_ORDER_WRITE_THRESHOLD).
    """
    cm5 = _parse_cm5_charges(lines)
    bo_mayer = _parse_mayer_bond_orders(lines)
    bo_total = _parse_total_bond_orders(lines)

    bo_raw = bo_mayer if bo_mayer else bo_total
    bo = {
        pair: v for pair, v in bo_raw.items()
        if v >= BOND_ORDER_WRITE_THRESHOLD
    }
    bo_source = "Mayer" if bo_mayer else ("Total" if bo_total else None)

    out = []
    if cm5:
        out.append(f"CM5_COUNT = {len(cm5)}")
        for idx in sorted(cm5):
            out.append(f"CM5 {idx} {cm5[idx]: .10f}")
    if bo:
        out.append(f"BO_SOURCE = {bo_source}")
        out.append(f"BO_COUNT = {len(bo)}")
        for (i, j), v in sorted(bo.items()):
            out.append(f"BO {i} {j} {v: .10f}")

    return out if out else None


def _parse_temperature_from_log(lines: List[str]):
    temps = []
    pattern = re.compile(r"^\s*T\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*K;", re.I)
    for line in lines:
        m = pattern.match(line)
        if not m:
            continue
        try:
            temps.append(float(m.group(1)))
        except Exception:
            continue
    if not temps:
        return None
    # Prefer the standard thermochemistry temperature when present.
    for t in temps:
        if abs(t - 298.15) < 1e-2:
            return t
    return temps[0]


def _parse_point_group_from_log(lines: List[str]):
    pg = None
    for line in lines:
        if "Full point group" in line:
            toks = line.split()
            if "group" in toks:
                try:
                    pg = toks[toks.index("group") + 1]
                except Exception:
                    pass
    return pg


def _build_rotational_meta_lines(lines: List[str]):
    out = []
    t = _parse_temperature_from_log(lines)
    if t is not None:
        out.append(f"T_K = {t:.2f}")
    pg = _parse_point_group_from_log(lines)
    if pg:
        out.append(f"group = {pg}")
    return out if out else None


def _build_rotational_lines_from_log(lines: List[str], deltabvib_lines=None):
    meta_lines = _build_rotational_meta_lines(lines) or []
    if not _has_output_pickett(lines):
        if deltabvib_lines:
            meta_lines.extend(deltabvib_lines)
        return meta_lines or None

    spfit = _parse_spfit_like_block(lines)
    rot_mhz = _parse_rot_constants_mhz(lines)
    asym = _parse_reduced_constants(lines, "Constants in the Asymmetrically reduced Hamiltonian")
    sym = _parse_reduced_constants(lines, "Constants in the Symmetrically Reduced Hamiltonian")

    lines_out = list(meta_lines)

    # Rotational constants (prefer SPFIT block if present, else A0/B0/C0, then A00/B00/C00, then Ae/Be/Ce)
    def _pick_rot(key_base):
        for k in (f"{key_base}0", f"{key_base}00", f"{key_base}e"):
            if k in rot_mhz:
                return rot_mhz[k]
        return None

    A = spfit.get("A")
    B = spfit.get("B")
    C = spfit.get("C")
    if A is None:
        A = _pick_rot("A")
    if B is None:
        B = _pick_rot("B")
    if C is None:
        C = _pick_rot("C")

    if A is not None:
        lines_out.append(f"A_MHz = {A:.13f}")
    if B is not None:
        lines_out.append(f"B_MHz = {B:.13f}")
    if C is not None:
        lines_out.append(f"C_MHz = {C:.13f}")

    # Distortion constants (S reduction, then A reduction)
    def _append_dist(keys, data):
        for k in keys:
            if k in data:
                lines_out.append(f"{k} = {data[k]:.12f}")

    # From SPFIT-like block (if present)
    spfit_map = {
        "DJ": "DJ_MHz",
        "DJK": "DJK_MHz",
        "DK": "DK_MHz",
        "d1": "d1_MHz",
        "d2": "d2_MHz",
        "DelJ": "DelJ_MHz",
        "DelJK": "DelJK_MHz",
        "DelK": "DelK_MHz",
        "delJ": "delJ_MHz",
        "delK": "delK_MHz",
    }
    spfit_norm = {}
    for k, v in spfit.items():
        if k in spfit_map:
            spfit_norm[spfit_map[k]] = v

    s_data = dict(sym)
    s_data.update(spfit_norm)
    a_data = dict(asym)
    a_data.update(spfit_norm)

    _append_dist(("DJ_MHz", "DJK_MHz", "DK_MHz", "d1_MHz", "d2_MHz"), s_data)
    _append_dist(("DelJ_MHz", "DelJK_MHz", "DelK_MHz", "delJ_MHz", "delK_MHz"), a_data)

    # Dipole components (Debye) from Pickett block, if available
    dip = _parse_pickett_dipole_debye(lines)
    if dip is not None:
        da, db, dc = dip
        lines_out.append(f"Dipole_a_D = {da:.8f}")
        lines_out.append(f"Dipole_b_D = {db:.8f}")
        lines_out.append(f"Dipole_c_D = {dc:.8f}")

    # Quadrupole coupling constants (Chi), if present
    chi = _parse_pickett_quadrupole_chi(lines)
    for k in sorted(chi.keys()):
        lines_out.append(f"{k} = {chi[k]:.10f}")

    # Sextic centrifugal distortion constants (Phi/H/phi/h), if present
    sextic = _parse_pickett_sextic_constants(lines)
    for k in sorted(sextic.keys()):
        lines_out.append(f"{k} = {sextic[k]:.10f}")

    # DeltaBvib inside #ROTATIONAL if available (MUST be last)
    if deltabvib_lines:
        lines_out.extend(deltabvib_lines)

    if not lines_out:
        return None
    return lines_out

def read_gaussian_properties(path):
    """
    Read Gaussian LOG/OUT and append properties (VIBRATIONAL/ROTATIONAL)
    into existing xyzin WITHOUT overwriting geometry.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() in {".fchk", ".fch"}:
        raise ValueError("Gaussian properties require .log/.out (fchk not supported here).")

    log_path, _ = _find_gaussian_companion_files(path)
    if log_path is None:
        raise ValueError(f"Unsupported Gaussian file: {path}")

    _copy_log_to_logfin(log_path)
    log_path = Path("working") / "logfin"

    lines = log_path.read_text(errors="ignore").splitlines()

    # --- DeltaBvib section (always present) ---
    vals, reason = _parse_deltabvib_from_log(log_path)
    if vals is None:
        deltabvib_lines = [f"NOT AVAILABLE: {reason}"]
        deltabvib_lines_for_rot = None
    else:
        dva, dvb, dvc = vals
        deltabvib_lines = [
            f"DVibA_MHz= {dva:.6f}",
            f"DVibB_MHz= {dvb:.6f}",
            f"DVibC_MHz= {dvc:.6f}",
        ]
        deltabvib_lines_for_rot = deltabvib_lines

    # --- VIBRATIONAL section from harmonic log data (if present) ---
    vib_lines = _build_vibrational_lines_from_log(lines)

    # --- ROTATIONAL section from output=pickett (if present) ---
    rot_lines = _build_rotational_lines_from_log(lines, deltabvib_lines=deltabvib_lines_for_rot)
    if not rot_lines and deltabvib_lines:
        rot_lines = deltabvib_lines
    gauss_topo_lines = _build_gaussian_topology_lines(lines)

    # --- Update xyzin sections only (do not touch geometry) ---
    if vib_lines:
        lines_xyzin = read_xyzin_lines()
        lines_xyzin = remove_section("VIBRATIONAL", lines_xyzin)
        write_xyzin_lines(lines_xyzin)
        append_section("VIBRATIONAL", vib_lines)

    if rot_lines:
        lines_xyzin = read_xyzin_lines()
        lines_xyzin = remove_section("DELTABVIB", lines_xyzin)
        lines_xyzin = remove_section("ROTATIONAL", lines_xyzin)
        write_xyzin_lines(lines_xyzin)
        append_section("ROTATIONAL", rot_lines)

    # --- GAUSSIAN_TOPOLOGY section (CM5 / bond orders) ---
    lines_xyzin = read_xyzin_lines()
    lines_xyzin = remove_section("GAUSSIAN_TOPOLOGY", lines_xyzin)
    write_xyzin_lines(lines_xyzin)
    if gauss_topo_lines:
        append_section("GAUSSIAN_TOPOLOGY", gauss_topo_lines)


def read_gaussian(path):
    """
    Read Gaussian output (.out/.log) or formatted checkpoint (.fchk/.fch)
    and update xyzin (XYZ + #BASIC + #DELTABVIB).

    Policy:
    - GUI provides a single file.
    - If a companion file with the same basename exists in the same directory,
      it is used automatically (.log/.out paired with .fchk/.fch).
    - Geometry and DeltaBvib are taken from the LOG when available.
    - If FCHK exists, it is copied to working/fchkin for downstream modules.
    - If LOG/OUT exists, it is copied to working/logfin and parsing happens there.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    log_path, fchk_path = _find_gaussian_companion_files(path)
    if log_path is None and fchk_path is None:
        raise ValueError(f"Unsupported Gaussian file: {path}")

    if log_path is None:
        print("WARNING: No Gaussian LOG/OUT file found (DeltaBvib will not be available).")
    if fchk_path is None:
        print("WARNING: No Gaussian FCHK/FCH file found (properties/normal modes may be unavailable).")

    # Copy companion files if available (then work ONLY on the copies)
    if fchk_path is not None:
        _copy_fchk_to_fchkin(fchk_path)
        fchk_path = Path("working") / "fchkin"
    if log_path is not None:
        _copy_log_to_logfin(log_path)
        log_path = Path("working") / "logfin"

    # --- Read geometry ---
    if log_path is not None:
        lines = log_path.read_text(errors="ignore").splitlines()
        charge, multiplicity = _parse_charge_mult(lines)
        geometries = _parse_orientations(lines)
        if not geometries:
            raise RuntimeError(f"No orientation found in Gaussian file: {log_path}")
        final_geom = geometries[-1]
        atoms = [a for a, _, _, _ in final_geom]
        coords = [(x, y, z) for _, x, y, z in final_geom]
        comment = log_path.name
    else:
        # FCHK-only mode
        _read_fchk_and_write_xyzin(fchk_path)
        lines_xyzin = read_xyzin_lines()
        lines_xyzin = remove_section("DELTABVIB", lines_xyzin)
        lines_xyzin = remove_section("VIBRATIONAL", lines_xyzin)
        lines_xyzin = remove_section("GAUSSIAN_TOPOLOGY", lines_xyzin)
        write_xyzin_lines(lines_xyzin)
        # Append DeltaBvib (not available) into #ROTATIONAL and return
        append_section("ROTATIONAL", ["NOT AVAILABLE: log file missing"])
        return

    # --- DeltaBvib section (always present) ---
    vals, reason = _parse_deltabvib_from_log(log_path)
    if vals is None:
        deltabvib_lines = [f"NOT AVAILABLE: {reason}"]
        deltabvib_lines_for_rot = None
    else:
        dva, dvb, dvc = vals
        deltabvib_lines = [
            f"DVibA_MHz= {dva:.6f}",
            f"DVibB_MHz= {dvb:.6f}",
            f"DVibC_MHz= {dvc:.6f}",
        ]
        deltabvib_lines_for_rot = deltabvib_lines

    _write_xyzin(atoms, coords, charge, multiplicity, comment)
    # Remove any legacy #DELTABVIB block (DeltaBvib now goes into #ROTATIONAL)
    lines_xyzin = read_xyzin_lines()
    lines_xyzin = remove_section("DELTABVIB", lines_xyzin)
    lines_xyzin = remove_section("VIBRATIONAL", lines_xyzin)
    write_xyzin_lines(lines_xyzin)

    # --- VIBRATIONAL section from harmonic log data (if present) ---
    vib_lines = _build_vibrational_lines_from_log(lines)
    if vib_lines:
        append_section("VIBRATIONAL", vib_lines)

    # --- ROTATIONAL section from output=pickett (if present) ---
    rot_lines = _build_rotational_lines_from_log(lines, deltabvib_lines=deltabvib_lines_for_rot)
    gauss_topo_lines = _build_gaussian_topology_lines(lines)
    if rot_lines:
        append_section("ROTATIONAL", rot_lines)
    elif deltabvib_lines:
        # Ensure DeltaBvib is not lost even without pickett
        append_section("ROTATIONAL", deltabvib_lines)

    # --- GAUSSIAN_TOPOLOGY section (CM5 / bond orders) ---
    lines_xyzin = read_xyzin_lines()
    lines_xyzin = remove_section("GAUSSIAN_TOPOLOGY", lines_xyzin)
    write_xyzin_lines(lines_xyzin)
    if gauss_topo_lines:
        append_section("GAUSSIAN_TOPOLOGY", gauss_topo_lines)


def _copy_fchk_to_fchkin(path: Path):
    """
    Copy the .fchk/.fch into the same folder where xyzin lives,
    using the canonical name 'fchkin'.
    """
    # xyzin_utils writes to working/xyzin -> that file is in parent folder of this module run context
    # safest is: put fchkin in current working dir (where Merlino runs)
    dst = Path("working") / "fchkin"
    dst.parent.mkdir(parents=True, exist_ok=True)

    try:
        shutil.copyfile(str(path), str(dst))
    except Exception:
        # fallback: try copy2 with metadata
        shutil.copy2(str(path), str(dst))


def _copy_log_to_logfin(path: Path):
    """
    Copy the .log/.out into the same folder where xyzin lives,
    using the canonical name 'logfin'.
    """
    dst = Path("working") / "logfin"
    dst.parent.mkdir(parents=True, exist_ok=True)

    try:
        shutil.copyfile(str(path), str(dst))
    except Exception:
        shutil.copy2(str(path), str(dst))


def _read_fchk_and_write_xyzin(path: Path):
    lines = path.read_text(errors="ignore").splitlines()

    atomic_numbers = []
    coords = []

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        if line.startswith("Atomic numbers"):
            count = int(line.split()[-1])
            i += 1
            while len(atomic_numbers) < count and i < n:
                atomic_numbers.extend(int(x) for x in lines[i].split())
                i += 1
            continue

        if line.startswith("Current cartesian coordinates"):
            count = int(line.split()[-1])
            i += 1
            while len(coords) < count and i < n:
                coords.extend(float(x) for x in lines[i].split())
                i += 1
            continue

        i += 1

    atoms = [_atomic_symbol(z) for z in atomic_numbers]

    xyz = []
    for i in range(0, len(coords), 3):
        xyz.append(
            (
                coords[i] * BOHR_TO_ANG,
                coords[i + 1] * BOHR_TO_ANG,
                coords[i + 2] * BOHR_TO_ANG,
            )
        )

    _write_xyzin(atoms, xyz, 0, 1, path.name)


# ==================================================
# WRITE XYzin
# ==================================================
def _write_xyzin(atoms, coords, charge, multiplicity, comment):
    """
    Write XYZ block and #BASIC section to xyzin (Merlino 3.0).
    """
    xyz_lines = [
        str(len(atoms)),
        comment,
    ]

    for a, (x, y, z) in zip(atoms, coords):
        xyz_lines.append(
            f"{a:2s} {x:15.8f} {y:15.8f} {z:15.8f}"
        )

    # --- XYZ
    replace_xyz_block(xyz_lines)

    # --- BASIC
    basic_lines = [
        f"charge {charge}",
        f"multiplicity {multiplicity}",
        "group c1",
        "T_K = 298.15",
        "P_ATM = 1.000000",
    ]
    append_section("BASIC", basic_lines)

    # NOTE: #DELTABVIB block is no longer written here; DeltaBvib is appended
    # at the end of #ROTATIONAL when available.
