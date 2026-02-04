# advanced/tests/test_provin_writer.py

from pathlib import Path
from advanced.provin_writer import ProvinWriter


def test_provin_writer(tmp_path):
    writer = ProvinWriter(
        workdir=tmp_path,
        title="Water molecule",
        charge=0,
        multiplicity=1,
        keywords=["GNIC", "G26", "SE"],
    )

    provin = writer.write()

    assert provin.exists()

    content = provin.read_text().splitlines()

    assert content[0] == "# GNIC G26 SE"
    assert content[1] == ""
    assert content[2] == "Water molecule"
    assert content[3] == ""
    assert content[4] == "0 1"
    assert len(content) == 5

