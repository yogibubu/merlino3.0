from __future__ import annotations

from pathlib import Path
import re


def extract_symmetry_section(report_text: str) -> str:
    """
    Extract the GLOBAL SYMMETRY section from topology.report text.
    Returns an empty string when no symmetry block is available.
    """
    if not report_text:
        return ""
    marker = "GLOBAL SYMMETRY"
    start = report_text.find(marker)
    if start < 0:
        return ""
    return report_text[start:].strip()


def load_symmetry_section(report_path: Path) -> str:
    if not report_path.exists():
        return ""
    try:
        text = report_path.read_text(encoding="utf-8")
    except Exception:
        return ""
    return extract_symmetry_section(text)


def parse_symmetry_overview(sym_text: str) -> dict:
    """
    Parse key symmetry metadata from the extracted symmetry section.
    """
    info = {
        "point_group": None,
        "atom_classes": [],
    }
    if not sym_text:
        return info

    for line in sym_text.splitlines():
        s = line.strip()
        if s.startswith("Point group:"):
            info["point_group"] = s.split(":", 1)[1].strip()
        elif s.startswith("class") and "[" in s and "]" in s:
            # Accept only atom-class lines before internal-parameter block.
            info["atom_classes"].append(s)
        elif s.startswith("EQUIVALENT INTERNAL-PARAMETER CLASSES"):
            break
    return info


def parse_equivalent_parameter_classes(sym_text: str) -> dict[str, list[tuple[str, str]]]:
    """
    Parse per-kind equivalent primitive classes from the symmetry section.
    Returns mapping: kind -> [(class_id, members), ...]
    """
    out: dict[str, list[tuple[str, str]]] = {}
    if not sym_text:
        return out

    in_block = False
    current_kind = None
    kind_header = re.compile(r"^(bond|angle|dihedral|out_of_plane|linear_bend):\s+\d+\s+equivalent classes$")
    class_line = re.compile(r"^class\s+(\d+):\s+(.*)$")

    for raw in sym_text.splitlines():
        s = raw.strip()
        if not in_block:
            if s.startswith("EQUIVALENT INTERNAL-PARAMETER CLASSES"):
                in_block = True
            continue
        if not s or set(s) <= {"-"}:
            continue
        m_kind = kind_header.match(s)
        if m_kind:
            current_kind = m_kind.group(1)
            out.setdefault(current_kind, [])
            continue
        m_cls = class_line.match(s)
        if m_cls and current_kind is not None:
            out[current_kind].append((m_cls.group(1), m_cls.group(2)))
    return out
