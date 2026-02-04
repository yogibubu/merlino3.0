from pathlib import Path
import shutil
import sys

from PySide6.QtWidgets import QApplication

# ------------------------------------------------------------
# Costanti di sessione
# ------------------------------------------------------------

ROOT = Path(__file__).resolve().parent
WORKING = ROOT / "working"
XYZIN = WORKING / "xyzin"

# ------------------------------------------------------------
# Bootstrap xyzin (MINIMALE, CANONICO)
# ------------------------------------------------------------

def write_initial_xyzin():
    WORKING.mkdir(exist_ok=True)

    content = """2
Merlino bootstrap
H   0.00000   0.00000   0.00000
H   0.00000   0.00000   0.74000

#BASIC
CHARGE              0
SPIN_MULTIPLICITY   1
POINT_GROUP         C1
REPRESENTATION      Ir
T_K =               298.15
P_ATM =             1.000000
"""
    XYZIN.write_text(content)

# ------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------

def cleanup_working(debug=False):
    """
    Clean the CONTENT of working/ without removing the directory itself.
    """
    if debug:
        return

    if not WORKING.exists():
        return

    for item in WORKING.iterdir():
        try:
            if item.is_file() or item.is_symlink():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
        except Exception:
            pass

# ------------------------------------------------------------
# Main entry
# ------------------------------------------------------------

def run(debug=False):
    # 1) bootstrap
    write_initial_xyzin()

    # 2) start GUI
    app = QApplication(sys.argv)

    from gui.main_window import MainWindow
    win = MainWindow(WORKING)
    win.show()

    exit_code = app.exec()

    # 3) cleanup
    cleanup_working(debug=debug)

    sys.exit(exit_code)

