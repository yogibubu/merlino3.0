"""
Interactive CLI to read a Gaussian file and run the full Merlino pipeline.
"""

from pathlib import Path
import sys
import math

from gui.gaussian import read_gaussian
from gui.project_manager import ProjectManager
from geometry.vib_anh import read_xyzin, direct_sum_dos, write_dos


def _prompt_path() -> Path:
    raw = input("Percorso file Gaussian (.log/.out/.fchk/.fch): ").strip()
    return Path(raw) if raw else None


def main(argv=None) -> int:
    argv = list(argv or [])
    root = Path(__file__).resolve().parent
    working = root / "working"
    working.mkdir(exist_ok=True)

    if argv:
        gpath = Path(argv[0]).expanduser()
    else:
        gpath = _prompt_path()

    if not gpath:
        print("Nessun file fornito.", file=sys.stderr)
        return 1

    if not gpath.exists():
        print(f"File non trovato: {gpath}", file=sys.stderr)
        return 1

    read_gaussian(gpath)

    pm = ProjectManager(working)
    pm.run_full_workflow()

    dos_vib = working / "dos_vib.dat"
    if not dos_vib.exists():
        try:
            vib = read_xyzin(str(working / "xyzin"))
            n = len(vib.omega_cm1)
            vmax = [6 for _ in range(n)]
            ncap = [10.0 for _ in range(n)]
            dos = direct_sum_dos(vib.omega_cm1, vib.chi_cm1, vmax, 8000.0, 50.0, ncap)
            dos_logg = {b: math.log(c) for b, c in dos.items() if c > 0}
            write_dos(str(dos_vib), dos_logg, 0.0, 50.0)
        except Exception:
            print("Pipeline completata: rotational, thermo, topology")
            print("Nota: DOS vibrazionale non generato, rovib non eseguito.")
            return 0

    pm.rovib_options = {"vib_dos": str(dos_vib)}
    pm.run("rovib")
    print("Pipeline completata: rotational, thermo, topology, rovib")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
