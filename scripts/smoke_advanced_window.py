from pathlib import Path
import sys
from PySide6.QtWidgets import QApplication

from advanced.advanced_window import AdvancedCalculationsWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)

    # directory di lavoro (deve contenere xyzin)
    workdir = Path("work")
    workdir.mkdir(exist_ok=True)

    # directory di progetto (nome usato come titolo)
    project_dir = Path("test_project")
    project_dir.mkdir(exist_ok=True)

    win = AdvancedCalculationsWindow(
        workdir=workdir,
        project_dir=project_dir
    )
    win.show()

    sys.exit(app.exec_())

