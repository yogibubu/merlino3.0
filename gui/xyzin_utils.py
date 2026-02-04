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

