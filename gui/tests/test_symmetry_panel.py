from gui.symmetry_panel import (
    extract_symmetry_section,
    parse_equivalent_parameter_classes,
    parse_symmetry_overview,
)


def test_extract_symmetry_section_returns_global_block():
    text = (
        "BOND ORDERS\n"
        "...\n"
        "\n"
        "GLOBAL SYMMETRY\n"
        "---------------\n"
        "Point group: Td\n"
        "Equivalent atom classes (1-based):\n"
        "  class  1: [2, 3, 4, 5]\n"
        "\n"
        "EQUIVALENT INTERNAL-PARAMETER CLASSES\n"
        "------------------------------------\n"
        "bond: 1 equivalent classes\n"
    )
    out = extract_symmetry_section(text)
    assert out.startswith("GLOBAL SYMMETRY")
    assert "Point group: Td" in out
    assert "EQUIVALENT INTERNAL-PARAMETER CLASSES" in out


def test_extract_symmetry_section_missing_marker():
    assert extract_symmetry_section("TOPOLOGY ONLY\n") == ""


def test_parse_symmetry_overview_and_classes():
    sym = (
        "GLOBAL SYMMETRY\n"
        "---------------\n"
        "Point group: Td\n"
        "Equivalent atom classes (1-based):\n"
        "  class  1: [2, 3, 4, 5]\n"
        "\n"
        "EQUIVALENT INTERNAL-PARAMETER CLASSES\n"
        "------------------------------------\n"
        "bond: 1 equivalent classes\n"
        "  class  1: ['(1-2)', '(1-3)', '(1-4)', '(1-5)']\n"
        "angle: 1 equivalent classes\n"
        "  class  1: ['(2-1-3)', '(2-1-4)']\n"
    )
    ov = parse_symmetry_overview(sym)
    assert ov["point_group"] == "Td"
    assert ov["atom_classes"]

    by_kind = parse_equivalent_parameter_classes(sym)
    assert "bond" in by_kind
    assert by_kind["bond"][0][0] == "1"
    assert "(1-2)" in by_kind["bond"][0][1]
