from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import shlex

from .xyzin_sections import has_section, read_sectioned_lines, replace_section, section_content


XYZIN_ISOTOPOLOGUES_SECTION = "ISOTOPOLOGUES"
XYZIN_ISOTOPOLOGUES_SCHEMA = "merlino.xyzin.isotopologues.v1"


@dataclass(frozen=True)
class XyzinIsotopologueRecord:
    label: str
    substitutions: dict[int, int] = field(default_factory=dict)
    rotational_MHz: tuple[float, float, float] | None = None
    deltavib_MHz: tuple[float, float, float] | None = None
    deltavib_source: str = "unspecified"
    deltavib_convention: str = "subtract"
    deltael_MHz: tuple[float, float, float] | None = None
    deltael_source: str = "unspecified"
    deltael_convention: str = "subtract"
    sigma_MHz: tuple[float, float, float] | None = None


def mass_number(value) -> int:
    text = str(value).strip()
    aliases = {"D": 2, "T": 3}
    if text.upper() in aliases:
        return aliases[text.upper()]
    return int(text)


def parse_substitutions(text: str) -> dict[int, int]:
    result: dict[int, int] = {}
    text = (text or "").strip()
    if not text:
        return result
    for chunk in text.split(";"):
        atom_text, isotope_text = chunk.split(":", 1)
        atom_index = int(atom_text.strip())
        isotope_a = mass_number(isotope_text.strip())
        if atom_index < 1:
            raise ValueError("Substitution atom indexes are one-based")
        result[atom_index] = isotope_a
    return result


def format_substitutions(substitutions: dict[int, int]) -> str:
    return ";".join(f"{atom}:{mass}" for atom, mass in sorted(substitutions.items()))


def sigma_text_from_weight(weight: float) -> str:
    if weight <= 0.0:
        return ""
    return f"{(1.0 / weight) ** 0.5:.12g}"


def has_xyzin_isotopologues(path: Path) -> bool:
    return has_section(Path(path), XYZIN_ISOTOPOLOGUES_SECTION)


def read_xyzin_isotopologue_records(path: Path) -> tuple[XyzinIsotopologueRecord, ...]:
    lines = section_content(read_sectioned_lines(Path(path)), XYZIN_ISOTOPOLOGUES_SECTION)
    return parse_xyzin_isotopologue_records(lines)


def write_xyzin_isotopologue_records(path: Path, records: tuple[XyzinIsotopologueRecord, ...]) -> Path:
    target = Path(path)
    replace_section(target, XYZIN_ISOTOPOLOGUES_SECTION, xyzin_isotopologue_section_lines(records))
    return target


def merge_xyzin_isotopologue_records(path: Path, records: tuple[XyzinIsotopologueRecord, ...]) -> Path:
    target = Path(path)
    existing = read_xyzin_isotopologue_records(target) if has_xyzin_isotopologues(target) else ()
    merged = _merge_records(existing, records)
    return write_xyzin_isotopologue_records(target, merged)


def parse_xyzin_isotopologue_records(lines: list[str]) -> tuple[XyzinIsotopologueRecord, ...]:
    records: list[XyzinIsotopologueRecord] = []
    current: dict[str, object] | None = None

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        tokens = shlex.split(line)
        if not tokens:
            continue
        keyword = tokens[0].upper()
        if keyword == "SCHEMA":
            schema = tokens[1] if len(tokens) > 1 else ""
            if schema != XYZIN_ISOTOPOLOGUES_SCHEMA:
                raise ValueError(f"Unsupported xyzin isotopologue schema: {schema!r}")
            continue
        if keyword in {"UNITS", "INDEXING"}:
            continue
        if keyword == "BEGIN":
            if current is not None:
                raise ValueError("Nested BEGIN in #ISOTOPOLOGUES")
            if len(tokens) < 2:
                raise ValueError("BEGIN in #ISOTOPOLOGUES needs a label")
            current = {
                "label": tokens[1],
                "substitutions": {},
                "rotational": None,
                "deltavib": None,
                "deltavib_source": "unspecified",
                "deltavib_convention": "subtract",
                "deltael": None,
                "deltael_source": "unspecified",
                "deltael_convention": "subtract",
                "sigma": None,
            }
            continue
        if keyword == "END":
            if current is None:
                raise ValueError("END without BEGIN in #ISOTOPOLOGUES")
            records.append(_record_from_mapping(current))
            current = None
            continue
        if current is None:
            raise ValueError(f"#ISOTOPOLOGUES entry outside BEGIN/END: {line}")
        if keyword == "DEFINITION":
            definition = " ".join(tokens[1:]).strip()
            current["substitutions"] = {} if not definition or definition.lower() == "parent" else parse_substitutions(definition)
            continue
        values = _assignment_dict(tokens[1:])
        if keyword == "ROTATIONAL_MHZ":
            current["rotational"] = _abc(values)
        elif keyword == "DELTAVIB_MHZ":
            current["deltavib"] = _abc_default(values)
            current["deltavib_source"] = str(values.get("SOURCE", "unspecified"))
            current["deltavib_convention"] = str(values.get("CONVENTION", "subtract"))
        elif keyword == "DELTAEL_MHZ":
            current["deltael"] = _abc_default(values)
            current["deltael_source"] = str(values.get("SOURCE", "unspecified"))
            current["deltael_convention"] = str(values.get("CONVENTION", "subtract"))
        elif keyword == "SIGMA_MHZ":
            sigmas = _abc(values)
            if any(sigma <= 0.0 for sigma in sigmas):
                raise ValueError("SIGMA_MHZ values in #ISOTOPOLOGUES must be positive")
            current["sigma"] = sigmas
        else:
            raise ValueError(f"Unsupported #ISOTOPOLOGUES keyword: {keyword}")

    if current is not None:
        raise ValueError("Unclosed BEGIN in #ISOTOPOLOGUES")
    if not records:
        raise ValueError("#ISOTOPOLOGUES contains no records")
    return tuple(records)


def xyzin_isotopologue_section_lines(records: tuple[XyzinIsotopologueRecord, ...]) -> list[str]:
    lines = [
        f"SCHEMA {XYZIN_ISOTOPOLOGUES_SCHEMA}",
        "UNITS ROTATIONAL=MHz DELTAVIB=MHz DELTAEL=MHz SIGMA=MHz",
        "INDEXING ATOMS=ONE_BASED",
    ]
    for record in records:
        label = record.label.strip() or "unnamed"
        definition = format_substitutions(record.substitutions) if record.substitutions else "parent"
        lines.extend(
            [
                f"BEGIN {shlex.quote(label)}",
                f"DEFINITION {definition}",
            ]
        )
        if record.rotational_MHz is not None:
            lines.append(_triple_line("ROTATIONAL_MHZ", record.rotational_MHz))
        if record.deltavib_MHz is not None:
            lines.append(
                _triple_line(
                    "DELTAVIB_MHZ",
                    record.deltavib_MHz,
                    source=record.deltavib_source,
                    convention=record.deltavib_convention,
                )
            )
        if record.deltael_MHz is not None:
            lines.append(
                _triple_line(
                    "DELTAEL_MHZ",
                    record.deltael_MHz,
                    source=record.deltael_source,
                    convention=record.deltael_convention,
                )
            )
        if record.sigma_MHz is not None:
            lines.append(_triple_line("SIGMA_MHZ", record.sigma_MHz))
        lines.append("END")
    return lines


def _merge_records(
    existing: tuple[XyzinIsotopologueRecord, ...],
    incoming: tuple[XyzinIsotopologueRecord, ...],
) -> tuple[XyzinIsotopologueRecord, ...]:
    ordered: list[XyzinIsotopologueRecord] = list(existing)
    by_label = {record.label: idx for idx, record in enumerate(ordered)}
    for record in incoming:
        idx = by_label.get(record.label)
        if idx is None:
            by_label[record.label] = len(ordered)
            ordered.append(record)
        else:
            ordered[idx] = _merge_record(ordered[idx], record)
    return tuple(ordered)


def _merge_record(old: XyzinIsotopologueRecord, new: XyzinIsotopologueRecord) -> XyzinIsotopologueRecord:
    return XyzinIsotopologueRecord(
        label=new.label or old.label,
        substitutions=new.substitutions,
        rotational_MHz=new.rotational_MHz if new.rotational_MHz is not None else old.rotational_MHz,
        deltavib_MHz=new.deltavib_MHz if new.deltavib_MHz is not None else old.deltavib_MHz,
        deltavib_source=new.deltavib_source if new.deltavib_MHz is not None else old.deltavib_source,
        deltavib_convention=new.deltavib_convention if new.deltavib_MHz is not None else old.deltavib_convention,
        deltael_MHz=new.deltael_MHz if new.deltael_MHz is not None else old.deltael_MHz,
        deltael_source=new.deltael_source if new.deltael_MHz is not None else old.deltael_source,
        deltael_convention=new.deltael_convention if new.deltael_MHz is not None else old.deltael_convention,
        sigma_MHz=new.sigma_MHz if new.sigma_MHz is not None else old.sigma_MHz,
    )


def _record_from_mapping(item: dict[str, object]) -> XyzinIsotopologueRecord:
    return XyzinIsotopologueRecord(
        label=str(item["label"]),
        substitutions=dict(item.get("substitutions", {})),
        rotational_MHz=item.get("rotational") if isinstance(item.get("rotational"), tuple) else None,
        deltavib_MHz=item.get("deltavib") if isinstance(item.get("deltavib"), tuple) else None,
        deltavib_source=str(item.get("deltavib_source", "unspecified")),
        deltavib_convention=str(item.get("deltavib_convention", "subtract")),
        deltael_MHz=item.get("deltael") if isinstance(item.get("deltael"), tuple) else None,
        deltael_source=str(item.get("deltael_source", "unspecified")),
        deltael_convention=str(item.get("deltael_convention", "subtract")),
        sigma_MHz=item.get("sigma") if isinstance(item.get("sigma"), tuple) else None,
    )


def _assignment_dict(tokens: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for token in tokens:
        if "=" not in token:
            raise ValueError(f"Expected KEY=VALUE assignment, got {token!r}")
        key, value = token.split("=", 1)
        values[key.strip().upper()] = value.strip()
    return values


def _abc(values: dict[str, str]) -> tuple[float, float, float]:
    return float(values["A"]), float(values["B"]), float(values["C"])


def _abc_default(values: dict[str, str]) -> tuple[float, float, float]:
    return float(values.get("A", 0.0)), float(values.get("B", 0.0)), float(values.get("C", 0.0))


def _triple_line(
    keyword: str,
    values: tuple[float, float, float],
    *,
    source: str | None = None,
    convention: str | None = None,
) -> str:
    line = f"{keyword} A={values[0]:.12g} B={values[1]:.12g} C={values[2]:.12g}"
    if source is not None:
        line += f" SOURCE={shlex.quote(source or 'unspecified')}"
    if convention is not None:
        line += f" CONVENTION={shlex.quote(convention or 'subtract')}"
    return line
