from __future__ import annotations

import math
import json
import hashlib
from pathlib import Path

from geometry.vib_anh import (
    read_xyzin,
    direct_sum_dos,
    write_dos,
    q_from_dos,
    read_rotational_block,
    _rotor_kind,
    _coerce_float,
    _coerce_int,
    rot_dos_logg,
    convolve_log_dos,
)


class DosController:
    def __init__(self, working_dir: Path):
        self.working_dir = working_dir
        self._cache_path = self.working_dir / "dos_cache.json"

    def run(self, xyzin_path: Path, settings: dict):
        if not xyzin_path.exists():
            return None

        vib_present, has_chi, neg_idx = self._vib_status(xyzin_path)
        if not vib_present:
            return {
                "input_type": "geometry-only",
                "vib_q": None,
                "rovib_q": None,
                "error": None,
            }

        is_ts = bool(neg_idx)
        ts_warn = None
        if is_ts and settings["dos_emin"] > 0.0:
            ts_warn = "TS detected: Emin should usually be 0 cm-1"

        if not settings["dos_do_vib"] and not settings["dos_do_rovib"]:
            return {
                "input_type": "ts" if is_ts else ("harmonic" if not has_chi else "anharmonic"),
                "vib_q": None,
                "rovib_q": None,
                "error": None,
            }

        cache_key = self._cache_key(xyzin_path, settings)
        cached = self._load_cache()
        if cached and cached.get("key") == cache_key:
            if self._cache_files_ok(settings):
                return cached.get("result")

        vib = read_xyzin(str(xyzin_path))
        omega, chi = self._filter_vib_data(vib, neg_idx)
        n = len(omega)
        vmax = [settings["dos_vmax"] for _ in range(n)]
        ncap = self._parse_float_list(settings["dos_ncap"], n)

        dos = direct_sum_dos(
            omega,
            chi,
            vmax,
            settings["dos_emax"],
            settings["dos_bin"],
            ncap,
        )

        dos_logg = {}
        for b, c in dos.items():
            if c <= 0:
                continue
            e_center = (b + 0.5) * settings["dos_bin"]
            if e_center < settings["dos_emin"] or e_center > settings["dos_emax"]:
                continue
            dos_logg[b] = math.log(c)

        if settings["dos_do_vib"]:
            dos_vib_path = self.working_dir / "dos_vib.dat"
            write_dos(str(dos_vib_path), dos_logg, 0.0, settings["dos_bin"])
            q_vib = self._q_from_binned(dos_logg, settings["dos_T"], settings["dos_bin"], 0.0)
            (self.working_dir / "vib_qt.dat").write_text(
                "# T_K Q_vib\n" + f"{settings['dos_T']:.6f} {q_vib:.12e}\n",
                encoding="utf-8",
            )
            if is_ts:
                n_vib = self._cumulative_logsumexp(dos_logg)
                write_dos(str(self.working_dir / "n_vib_ts.dat"), n_vib, 0.0, settings["dos_bin"])
        else:
            q_vib = None

        if settings["dos_do_rovib"]:
            rot_block = read_rotational_block(str(xyzin_path))
            rotor_kind = _rotor_kind(rot_block.get("rotor_type"))
            if rotor_kind is None:
                return {
                    "input_type": "ts" if is_ts else ("anharmonic" if has_chi else "harmonic"),
                    "vib_q": q_vib,
                    "rovib_q": None,
                    "error": self._merge_warning(ts_warn, "rovib skipped: missing rotor_type"),
                }
            A = _coerce_float(rot_block.get("a_mhz"))
            B = _coerce_float(rot_block.get("b_mhz"))
            C = _coerce_float(rot_block.get("c_mhz"))
            if A is None or B is None or C is None:
                return {
                    "input_type": "ts" if is_ts else ("anharmonic" if has_chi else "harmonic"),
                    "vib_q": q_vib,
                    "rovib_q": None,
                    "error": self._merge_warning(ts_warn, "rovib skipped: missing rotational constants"),
                }
            sigma = _coerce_int(rot_block.get("symm. number"))
            if sigma is None:
                sigma = _coerce_int(rot_block.get("sigma")) or 1

            emax_rot = settings["dos_emax_rot"] if settings["dos_emax_rot"] else settings["dos_emax"]
            jmax = settings["dos_jmax"] if settings["dos_jmax"] else None

            rot_logg = rot_dos_logg(
                rotor_kind, A, B, C, sigma, emax_rot, settings["dos_bin"], jmax=jmax
            )
            rovib_logg = convolve_log_dos(dos_logg, rot_logg)
            dos_rovib_path = self.working_dir / "dos_rovib.dat"
            write_dos(str(dos_rovib_path), rovib_logg, 0.0, settings["dos_bin"])
            q_rovib = self._q_from_binned(rovib_logg, settings["dos_T"], settings["dos_bin"], 0.0)
            (self.working_dir / "rovib_qt.dat").write_text(
                "# T_K Q_rovib\n" + f"{settings['dos_T']:.6f} {q_rovib:.12e}\n",
                encoding="utf-8",
            )
            if is_ts:
                n_rovib = self._cumulative_logsumexp(rovib_logg)
                write_dos(
                    str(self.working_dir / "n_rovib_ts.dat"), n_rovib, 0.0, settings["dos_bin"]
                )
        else:
            q_rovib = None

        result = {
            "input_type": "ts" if is_ts else ("anharmonic" if has_chi else "harmonic"),
            "vib_q": q_vib,
            "rovib_q": q_rovib,
            "error": ts_warn,
        }
        self._save_cache(cache_key, result)
        return result

    def _cache_key(self, xyzin_path: Path, settings: dict) -> str:
        try:
            content = xyzin_path.read_bytes()
        except Exception:
            content = b""
        settings_blob = json.dumps(settings, sort_keys=True, default=str).encode("utf-8")
        h = hashlib.sha256()
        h.update(content)
        h.update(b"\0")
        h.update(settings_blob)
        return h.hexdigest()

    def _cache_files_ok(self, settings: dict) -> bool:
        required = [self.working_dir / "dos_vib.dat", self.working_dir / "vib_qt.dat"]
        if settings.get("dos_do_rovib"):
            required.extend(
                [
                    self.working_dir / "dos_rovib.dat",
                    self.working_dir / "rovib_qt.dat",
                ]
            )
        for path in required:
            if not path.exists():
                return False
        return True

    def _load_cache(self):
        if not self._cache_path.exists():
            return None
        try:
            return json.loads(self._cache_path.read_text())
        except Exception:
            return None

    def _save_cache(self, key: str, result: dict):
        data = {"key": key, "result": result}
        try:
            self._cache_path.write_text(json.dumps(data, indent=2))
        except Exception:
            pass

    def _vib_status(self, xyzin_path: Path):
        lines = xyzin_path.read_text().splitlines()
        in_vib = False
        has_vib = False
        has_chi = False
        freqs = []
        for line in lines:
            s = line.strip()
            if s.startswith("#"):
                in_vib = (s.upper() == "#VIBRATIONAL")
                continue
            if not in_vib:
                continue
            if not s:
                continue
            has_vib = True
            if s.lower().startswith("chi_cm1"):
                has_chi = True
            if s.lower().startswith("freq") or s.lower().startswith("frequencies"):
                freqs += self._grab_numbers(s)
            elif freqs and not s[0].isalpha():
                freqs += self._grab_numbers(s)
        neg_idx = [i for i, f in enumerate(freqs) if f < 0.0]
        return has_vib, has_chi, neg_idx

    def _grab_numbers(self, s: str):
        import re
        return [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]

    def _filter_vib_data(self, vib, neg_idx):
        if not neg_idx:
            return vib.omega_cm1, vib.chi_cm1
        keep = [i for i in range(len(vib.omega_cm1)) if i not in set(neg_idx)]
        omega = [vib.omega_cm1[i] for i in keep]
        chi = [[vib.chi_cm1[i][j] for j in keep] for i in keep]
        return omega, chi

    def _cumulative_logsumexp(self, dos_logg):
        bins = sorted(dos_logg.keys())
        out = {}
        acc = None
        for b in bins:
            acc = self._logsumexp(acc, dos_logg[b])
            out[b] = acc
        return out

    def _logsumexp(self, a, b):
        if a is None:
            return b
        if b is None:
            return a
        m = a if a > b else b
        return m + math.log(math.exp(a - m) + math.exp(b - m))

    def _merge_warning(self, base, extra):
        if base and extra:
            return base + "; " + extra
        return base or extra

    def _parse_float_list(self, text: str, n: int):
        if not text:
            return None
        parts = [p for p in text.split(",") if p.strip()]
        if not parts:
            return None
        if len(parts) == 1:
            return [float(parts[0]) for _ in range(n)]
        if len(parts) != n:
            raise ValueError("ncap list length must match number of modes")
        return [float(x) for x in parts]

    def _q_from_binned(self, dos_logg, t_k, bin_cm1, emin_cm1):
        dos_e = {}
        for b, lg in dos_logg.items():
            e = emin_cm1 + (b + 0.5) * bin_cm1
            dos_e[e] = lg
        return q_from_dos(dos_e, t_k)
