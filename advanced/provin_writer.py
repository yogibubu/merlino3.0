from pathlib import Path
from typing import List


class ProvinWriter:
    """
    Writer for the 'provin' input file.

    Format (strict):

    # <keywords on a single line>

    <title>

    <charge> <multiplicity>
    """

    def __init__(
        self,
        workdir: Path,
        title: str,
        charge: int,
        multiplicity: int,
        keywords: List[str],
    ):
        self.workdir = Path(workdir)
        self.title = title.strip()
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        self.keywords = [kw.strip() for kw in keywords if kw.strip()]

    # ------------------------------------------------------------------

    def write(self) -> Path:
        """
        Write the provin file and return its path.
        """
        # Ensure working directory exists
        self.workdir.mkdir(parents=True, exist_ok=True)

        provin_path = self.workdir / "provin"

        with provin_path.open("w") as f:
            # line 1: keywords
            f.write("# " + " ".join(self.keywords) + "\n")

            # line 2: blank
            f.write("\n")

            # line 3: title
            f.write(self.title + "\n")

            # line 4: blank
            f.write("\n")

            # line 5: charge and multiplicity
            f.write(f"{self.charge} {self.multiplicity}\n")

        return provin_path

