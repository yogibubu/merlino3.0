from pathlib import Path

# --------------------------------------------------
# Paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
WORKING = ROOT / "working"
XYZIN = WORKING / "xyzin"


# --------------------------------------------------
# Low-level I/O
# --------------------------------------------------

def read_xyzin_lines():
    if not XYZIN.exists():
        return []
    return XYZIN.read_text().splitlines()


def write_xyzin_lines(lines):
    XYZIN.write_text("\n".join(lines) + "\n")


# --------------------------------------------------
# Section utilities
# --------------------------------------------------

def remove_section(section_name, lines):
    """
    Remove a section (#SECTION_NAME) and its content.
    """
    section_name = section_name.upper()
    out = []
    skip = False

    for line in lines:
        if line.strip().upper() == f"#{section_name}":
            skip = True
            continue

        if skip:
            if line.startswith("#"):
                skip = False
                out.append(line)
            continue

        out.append(line)

    return out


def append_section(section_name, content_lines):
    """
    Append a section at the end of xyzin, ensuring
    a single blank line before the section header.
    """
    lines = read_xyzin_lines()
    lines = remove_section(section_name, lines)

    # Ensure exactly one blank line before new section
    if lines and lines[-1].strip():
        lines.append("")

    lines.append(f"#{section_name.upper()}")
    for line in content_lines:
        lines.append(line)

    write_xyzin_lines(lines)


def set_rotational_keyvals(keyvals):
    """
    Ensure #ROTATIONAL exists and set/update key=value lines, preserving
    unrelated entries. Unknown keys are allowed on purpose so external
    rovibrational tools can attach compatibility metadata without breaking
    Merlino readers.
    """
    lines = read_xyzin_lines()
    if not lines:
        return False

    normalized = [(str(k).strip(), str(v).strip()) for k, v in keyvals.items()]
    if not normalized:
        return True

    keys_upper = {k.upper() for k, _ in normalized}
    new_lines = [f"{k}= {v}" for k, v in normalized]

    rot_start = None
    for i, line in enumerate(lines):
        if line.strip().upper() == "#ROTATIONAL":
            rot_start = i
            break

    if rot_start is None:
        append_section("ROTATIONAL", new_lines)
        return True

    out = []
    in_rot = False
    inserted = False
    for line in lines:
        if line.strip().upper() == "#ROTATIONAL":
            in_rot = True
            out.append(line)
            continue
        if in_rot:
            if line.startswith("#"):
                in_rot = False
                if not inserted:
                    out.extend(new_lines)
                    inserted = True
                out.append(line)
                continue

            stripped = line.strip()
            if "=" in stripped:
                k = stripped.split("=", 1)[0].strip().upper()
                if k in keys_upper:
                    continue
        out.append(line)

    if in_rot and not inserted:
        out.extend(new_lines)

    write_xyzin_lines(out)
    return True


def set_dvib_in_rotational(dva, dvb, dvc):
    """
    Ensure #ROTATIONAL exists and set DVibA/B/C lines, preserving other entries.
    """
    return set_rotational_keyvals(
        {
            "DVibA_MHz": f"{dva:.6f}",
            "DVibB_MHz": f"{dvb:.6f}",
            "DVibC_MHz": f"{dvc:.6f}",
        }
    )


# --------------------------------------------------
# XYZ block replacement
# --------------------------------------------------

def replace_xyz_block(xyz_lines):
    """
    Replace the XYZ block at the beginning of xyzin.
    """
    lines = read_xyzin_lines()

    # Skip existing XYZ block
    i = 0
    if lines:
        try:
            n = int(lines[0])
            i = n + 2
        except Exception:
            i = 0

    new_lines = []
    new_lines.extend(xyz_lines)

    # Preserve rest of file
    new_lines.extend(lines[i:])

    write_xyzin_lines(new_lines)
