from __future__ import annotations

from pathlib import Path


class XyzinService:
    def _remove_smiles_from_lines(self, lines):
        out, skip = [], False
        for line in lines:
            if line.strip().upper() == "#SMILES":
                skip = True
                continue
            if skip and not line.startswith("#"):
                continue
            skip = False
            out.append(line)
        return out

    def remove_smiles_section(self, xyzin_path: Path):
        if not xyzin_path.exists():
            return
        lines = xyzin_path.read_text().splitlines()
        out = self._remove_smiles_from_lines(lines)
        xyzin_path.write_text("\n".join(out) + "\n")

    def parse_basic_section(self, xyzin_path: Path):
        lines = xyzin_path.read_text().splitlines()
        if "#BASIC" not in lines:
            return {}
        i = lines.index("#BASIC") + 1
        basic = {}
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith("#"):
                break
            toks = line.split()
            if len(toks) >= 2:
                basic[toks[0].upper()] = toks[-1]
            i += 1
        return basic

    def format_basic_line(self, key: str, basic_values: dict):
        if key == "CHARGE":
            return f"CHARGE              {basic_values['basic_charge']}"
        if key == "SPIN_MULTIPLICITY":
            return f"SPIN_MULTIPLICITY   {basic_values['basic_mult']}"
        if key == "POINT_GROUP":
            return f"POINT_GROUP         {basic_values['basic_pg']}"
        if key == "REPRESENTATION":
            return f"REPRESENTATION      {basic_values['basic_representation']}"
        if key == "T_K":
            return f"T_K =               {basic_values['basic_T']:.2f}"
        if key == "P_ATM":
            return f"P_ATM =             {basic_values['basic_P_atm']:.6f}"
        return ""

    def _update_basic_lines(self, lines, gui_keys, basic_values: dict):
        try:
            i0 = next(i for i, ln in enumerate(lines) if ln.strip().upper() == "#BASIC")
        except StopIteration:
            lines.append("#BASIC")
            i0 = len(lines) - 1

        i = i0 + 1
        while i < len(lines) and lines[i].strip() and not lines[i].startswith("#"):
            i += 1

        new_block = []
        for key in gui_keys:
            new_block.append(self.format_basic_line(key, basic_values))

        return lines[:i0 + 1] + new_block + lines[i:]

    def update_basic_section(self, xyzin_path: Path, gui_keys, basic_values: dict):
        if not xyzin_path.exists():
            return
        lines = xyzin_path.read_text().splitlines()
        out = self._update_basic_lines(lines, gui_keys, basic_values)
        xyzin_path.write_text("\n".join(out) + "\n")

    def update_post_input(self, xyzin_path: Path, gui_keys, basic_values: dict, remove_smiles: bool):
        if not xyzin_path.exists():
            return
        lines = xyzin_path.read_text().splitlines()
        if remove_smiles:
            lines = self._remove_smiles_from_lines(lines)
        lines = self._update_basic_lines(lines, gui_keys, basic_values)
        xyzin_path.write_text("\n".join(lines) + "\n")

        try:
            i0 = next(i for i, ln in enumerate(lines) if ln.strip().upper() == "#BASIC")
        except StopIteration:
            lines.append("#BASIC")
            i0 = len(lines) - 1

        i = i0 + 1
        while i < len(lines) and lines[i].strip() and not lines[i].startswith("#"):
            i += 1

        new_block = []
        for key in gui_keys:
            new_block.append(self.format_basic_line(key, basic_values))

        out = lines[:i0 + 1] + new_block + lines[i:]
        xyzin_path.write_text("\n".join(out) + "\n")

    def summary(self, xyzin_path: Path):
        try:
            basic = self.parse_basic_section(xyzin_path)
            return (
                f"Charge: {basic.get('CHARGE', '?')}\n"
                f"Multiplicity: {basic.get('SPIN_MULTIPLICITY', '?')}\n"
                f"Point group: {basic.get('POINT_GROUP', '?')}\n"
                f"Representation: {basic.get('REPRESENTATION', 'Ir')}\n"
                f"T (K): {basic.get('T_K', '?')}\n"
                f"P (atm): {basic.get('P_ATM', '?')}"
            )
        except Exception:
            return "BASIC section unreadable."

    def has_rotational_block(self, xyzin_path: Path) -> bool:
        if not xyzin_path.exists():
            return False
        lines = xyzin_path.read_text().splitlines()
        return any(line.strip().upper() == "#ROTATIONAL" for line in lines)
