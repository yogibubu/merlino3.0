"""
Pytest module for Gaussian reader (Merlino 3.0).

Test files are located in:
    gui/tests/gaussian/
"""

from pathlib import Path

from gui.gaussian import read_gaussian


# --------------------------------------------------
# Paths (relative to this file)
# --------------------------------------------------
HERE = Path(__file__).resolve().parent
GAUSSIAN_TEST_DIR = HERE / "tests" / "gaussian"


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def init_xyzin(working_dir: Path):
    xyzin = working_dir / "xyzin"
    xyzin.write_text(
        "2\n"
        "H2\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 0.7\n\n"
        "#BASIC\n"
        "charge 0\n"
        "multiplicity 1\n"
        "group c1\n"
    )
    return xyzin


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


def check_xyz(xyzin_text: str):
    lines = xyzin_text.splitlines()
    nat = int(lines[0])
    assert nat > 0
    assert len(lines) >= nat + 2


# --------------------------------------------------
# TEST
# --------------------------------------------------
def test_gaussian_outputs(tmp_path):
    assert GAUSSIAN_TEST_DIR.exists(), (
        f"Gaussian test directory not found: {GAUSSIAN_TEST_DIR}"
    )

    gaussian_files = [
        f for f in GAUSSIAN_TEST_DIR.iterdir()
        if f.suffix.lower() in {".out", ".log", ".fchk", ".fch"}
    ]

    assert gaussian_files, "No Gaussian test files found"

    for gfile in gaussian_files:
        working = tmp_path / gfile.stem
        working.mkdir(exist_ok=True)

        xyzin = init_xyzin(working)

        # Run reader
        read_gaussian(gfile)

        # Validate xyzin
        text = xyzin.read_text()
        check_xyz(text)

        basic = parse_basic(text)
        assert "charge" in basic
        assert "multiplicity" in basic

        print(
            f"OK {gfile.name}: "
            f"charge={basic['charge']} "
            f"multiplicity={basic['multiplicity']}"
        )

