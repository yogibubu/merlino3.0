"""
Pytest module for Z-matrix reader (Merlino 3.0).

Test files are located in:
    gui/tests/zmat/
"""

from pathlib import Path

from gui.zmat_reader import read_zmat


# --------------------------------------------------
# Paths (relative to this file)
# --------------------------------------------------
HERE = Path(__file__).resolve().parent
ZMAT_TEST_DIR = HERE / "tests" / "zmat"
ROOT = HERE.parent


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def init_xyzin(working_dir: Path):
    """
    Initialize a minimal valid xyzin file.
    """
    xyzin = working_dir / "xyzin"
    xyzin.write_text(
        "1\n"
        "Dummy\n"
        "He 0.0 0.0 0.0\n\n"
        "#BASIC\n"
        "charge 0\n"
        "multiplicity 1\n"
        "group c1\n"
    )
    return xyzin


def _snapshot(path: Path):
    return path.read_bytes() if path.exists() else None


def _restore(path: Path, payload):
    if payload is None:
        if path.exists():
            path.unlink()
    else:
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(payload)


def check_xyz(xyzin_text: str):
    """
    Minimal sanity check for XYZ block.
    """
    lines = xyzin_text.splitlines()
    nat = int(lines[0])
    assert nat > 0
    assert len(lines) >= nat + 2


def parse_basic(xyzin_text: str):
    """
    Minimal local parser for #BASIC section.
    """
    lines = xyzin_text.splitlines()
    if "#BASIC" not in lines:
        return {}

    i = lines.index("#BASIC") + 1
    basic = {}

    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            break
        key, value = line.split()
        basic[key] = value
        i += 1

    return basic


# --------------------------------------------------
# TEST
# --------------------------------------------------
def test_zmat_outputs(tmp_path):
    assert ZMAT_TEST_DIR.exists(), (
        f"Z-matrix test directory not found: {ZMAT_TEST_DIR}"
    )

    zmat_files = [
        f for f in ZMAT_TEST_DIR.iterdir()
        if f.suffix.lower() in {".zmat", ".z", ".inp"}
    ]

    assert zmat_files, "No Z-matrix test files found"

    global_xyzin = ROOT / "working" / "xyzin"
    global_xyzin_snapshot = _snapshot(global_xyzin)
    try:
        for zfile in zmat_files:
            working = tmp_path / zfile.stem
            working.mkdir(exist_ok=True)

            xyzin = init_xyzin(working)

            # Run reader
            read_zmat(zfile)

            # Validate xyzin
            text = xyzin.read_text()

            # XYZ block must exist and be valid
            check_xyz(text)

            # BASIC section must be preserved
            basic = parse_basic(text)
            assert "charge" in basic
            assert "multiplicity" in basic

            print(
                f"OK {zfile.name}: "
                f"nat={text.splitlines()[0]} "
                f"charge={basic['charge']} "
                f"multiplicity={basic['multiplicity']}"
            )
    finally:
        _restore(global_xyzin, global_xyzin_snapshot)
