from pathlib import Path

from gui.readers import read_structure


class DummyInput:
    def __init__(self, kind, **kwargs):
        self.kind = kind
        for k, v in kwargs.items():
            setattr(self, k, v)


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def get_project_root():
    """
    Project root = parent of gui/
    """
    return Path(__file__).resolve().parents[1]


def init_working():
    """
    Initialize ROOT/working/xyzin according to Merlino contract.
    """
    root = get_project_root()
    working = root / "working"
    working.mkdir(exist_ok=True)

    xyzin = working / "xyzin"
    xyzin.write_text(
        "2\n"
        "H2\n"
        "H 0 0 0\n"
        "H 0 0 0.7\n\n"
        "#BASIC\n"
        "charge 0\n"
        "multiplicity 1\n"
        "group c1\n"
    )
    return xyzin


# --------------------------------------------------
# TEST
# --------------------------------------------------
def test_xyz_reader():
    xyzin = init_working()

    xyz = get_project_root() / "test.xyz"
    xyz.write_text(
        "3\n"
        "water\n"
        "O 0 0 0\n"
        "H 0 0 1\n"
        "H 1 0 0\n"
    )

    try:
        inp = DummyInput(kind="xyz", path=xyz)
        read_structure(inp)

        content = xyzin.read_text()
        assert content.startswith("3\n")
        assert "O" in content

    finally:
        # cleanup
        if xyz.exists():
            xyz.unlink()

