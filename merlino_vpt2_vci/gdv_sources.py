from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


DEFAULT_GDV_SOURCE_ROOT = Path.home() / "gdv_j32p" / "gdv"


@dataclass(frozen=True)
class FortranDeck:
    """Location of a named `*Deck` section inside a GDV Fortran source."""

    name: str
    source: Path
    line: int


@dataclass(frozen=True)
class GDVVPT2VCISources:
    """Relevant GDV source files for the Merlino4 VPT2/VCI refactor."""

    root: Path
    vci_driver: Path | None
    anharmonic_data: Path | None
    utility_support: Path | None
    decks: tuple[FortranDeck, ...]


VCI_DECKS = (
    "VCI1MI",
    "VCIDrv",
    "VCIGen",
    "VCIInt",
    "VCIPrt",
    "VCIPT2",
    "VCIVar",
    "VPTCrs",
    "VPTRHS",
    "VPTWFC",
)

ANHARMONIC_DECKS = (
    "AnhFIO",
    "DNICFq",
    "DNICGF",
    "VibAlpha",
    "VibCNM",
    "VibCor",
    "VPT2En",
)

UTILITY_DECKS = (
    "NHDiag",
    "VibFq2",
    "VibFrq",
    "VibOvl",
    "VibPrj",
    "VibSym",
    "VibTbl",
)


def _scan_decks(source: Path, wanted: set[str]) -> tuple[FortranDeck, ...]:
    if not source.exists():
        return ()
    decks: list[FortranDeck] = []
    for line_no, line in enumerate(source.read_text(errors="ignore").splitlines(), 1):
        if not line.startswith("*Deck"):
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[1] in wanted:
            decks.append(FortranDeck(parts[1], source, line_no))
    return tuple(decks)


def discover_gdv_vpt2_vci_sources(
    root: Path | None = None,
) -> GDVVPT2VCISources:
    """Find the GDV source files that contain VPT2/VCI and support routines."""
    gdv_root = Path(root) if root is not None else DEFAULT_GDV_SOURCE_ROOT
    vci_driver = gdv_root / "l717.F"
    anharmonic_data = gdv_root / "dinautil.F"
    utility_support = gdv_root / "utilnz.F"

    decks = (
        *_scan_decks(vci_driver, set(VCI_DECKS)),
        *_scan_decks(anharmonic_data, set(ANHARMONIC_DECKS)),
        *_scan_decks(utility_support, set(UTILITY_DECKS)),
    )

    return GDVVPT2VCISources(
        root=gdv_root,
        vci_driver=vci_driver if vci_driver.exists() else None,
        anharmonic_data=anharmonic_data if anharmonic_data.exists() else None,
        utility_support=utility_support if utility_support.exists() else None,
        decks=decks,
    )

