from pathlib import Path
import numpy as np

from survibfit.vibrational_internal import modes_from_gaussian_log


def test_vib_from_gaussian_log(tmp_path):
    log_path = Path("/Users/vincenzobarone/merlino3.0/gui/tests/gaussian/h2o.out")
    if not log_path.exists():
        return

    # Minimal fchk with atomic masses
    fchk = tmp_path / "h2o.fchk"
    fchk.write_text(
        "Atomic masses                            3\n"
        "  15.994915 1.007825 1.007825\n"
    )

    freqs, modes_q, U, prims = modes_from_gaussian_log(log_path, fchk_path=fchk)
    assert len(freqs) > 0
    assert modes_q.shape[0] == len(freqs)
    assert U.shape[0] == len(prims)
