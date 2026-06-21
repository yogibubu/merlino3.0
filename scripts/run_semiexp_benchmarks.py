from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from merlino_semiexp import SemiexperimentalFitRequest, read_observations, read_semiexperimental_job
from merlino_semiexp.fit import _atomic_number, _gic_model, _gicforge_a1_mask, _make_gicforge_backend
from merlino_semiexp.geometry_input import read_geometry_input
from merlino_semiexp.report import SemiexperimentalBenchmarkCase, benchmark_csv, run_semiexperimental_benchmark


GEOMETRIES = {
    "pyridine": """11
pyridine benchmark, C2v topology
H     -2.494015    0.000000    0.000000
C     -1.410500    0.000000    0.000000
C     -0.699065    1.194769    0.000000
H     -1.206310    2.151475    0.000000
C      0.692105    1.140055    0.000000
H      1.275926    2.054914    0.000000
N      1.390142    0.000000    0.000000
C      0.692105   -1.140055    0.000000
H      1.275926   -2.054914    0.000000
C     -0.699065   -1.194769    0.000000
H     -1.206310   -2.151475    0.000000
""",
    "uracil": """12
uracil benchmark
O      2.762610    0.073989    0.143605
C      1.382866    0.218650    0.054242
N      0.841722    1.456294   -0.095514
C     -0.503598    1.623167   -0.185155
C     -1.348077    0.515041   -0.123611
C     -0.776884   -0.748605    0.030311
O     -1.589209   -1.875729    0.095457
N      0.574465   -0.873514    0.116205
H      3.182936   -0.840906    0.255416
H     -0.916835    2.615835   -0.304133
H     -2.421614    0.635131   -0.193873
H     -1.188382   -2.799353    0.207051
""",
    "guanine": """16
guanine benchmark
N    -3.097661   -0.959464    0.0000
C    -1.752117   -0.483039   -0.0000
N    -0.718102   -1.380243   -0.0000
C     0.546223   -0.896522   -0.0000
N     1.682616   -1.598262   -0.0000
C     2.643469   -0.639746    0.0000
N     2.150071    0.630065    0.0000
C     0.820786    0.425058    0.0000
C    -0.215097    1.347250    0.0000
O     0.014202    2.586108   -0.0000
N    -1.497545    0.869919   -0.0000
H    -3.900146   -0.291510    0.0000
H    -3.293479   -1.985037    0.0000
H     1.804423   -2.635132    0.0000
H     3.702363   -0.862162   -0.0000
H    -2.295099    1.545612   -0.0000
""",
    "cyclopentadiene": """11
cyclopentadiene benchmark
C    6.40606548E-08    1.21099405E+00   -4.37769374E-05
C    1.16966652E+00    3.42104930E-01    3.60217132E-05
C   -1.16966649E+00    3.42105050E-01    3.19263268E-05
C    7.31740874E-01   -9.75616979E-01   -1.77891595E-05
C   -7.31740976E-01   -9.75616904E-01   -2.03514674E-05
H    2.26478781E+00    3.74863776E-01    1.27414604E-04
H   -2.26478777E+00    3.74864010E-01    1.20085357E-04
H    1.26945754E+00   -1.92001108E+00    6.26767155E-06
H   -1.26945774E+00   -1.92001095E+00    1.85908135E-06
H    5.85394489E-05    1.87871494E+00   -8.68912000E-01
H   -5.83447730E-05    1.87871720E+00    8.68822706E-01
""",
}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Merlino semiexperimental benchmark cases.")
    parser.add_argument("--outdir", type=Path, default=Path("benchmarks/semiexp_runs"))
    parser.add_argument("--max-iter", type=int, default=None)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    repo = Path(__file__).resolve().parents[1]
    geometry_paths = {
        "water": repo / "examples/semiexp/water/parent.xyz",
        "azulene": repo / "working/semiexp/azulene/parent.xyz",
    }
    for name, text in GEOMETRIES.items():
        path = outdir / f"{name}.xyz"
        path.write_text(text, encoding="utf-8")
        geometry_paths[name] = path

    gic_rows = []
    for name, path in sorted(geometry_paths.items()):
        try:
            geometry = read_geometry_input(path)
            atoms = tuple(geometry.atoms)
            coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
            z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
            backend = _make_gicforge_backend(atoms, outdir / "gicforge" / name)
            _prims, _u_matrix, labels = _gic_model(coords, z_numbers, backend=backend)
            a1_count = int(np.sum(_gicforge_a1_mask(labels)))
            gic_rows.append((name, "ok", len(atoms), backend.point_group or "UNKNOWN", len(labels), a1_count, backend.counter, ""))
        except Exception as exc:
            gic_rows.append((name, "error", 0, "UNKNOWN", 0, 0, 0, str(exc).replace("\n", " ")))

    with (outdir / "semiexp_gic_benchmarks.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "status", "natoms", "point_group", "n_gics", "n_a1_gics", "gicforge_calls", "message"])
        writer.writerows(gic_rows)

    azulene_job = read_semiexperimental_job(repo / "working/semiexp/azulene/azulene_msr_like.mse.toml")
    azulene_request = SemiexperimentalFitRequest(
        azulene_job.path,
        read_observations(azulene_job.observations),
        fixed_parameters=azulene_job.fixed_parameters,
        observable=azulene_job.observable,
        rotational_components=azulene_job.rotational_components,
        qm_predicates=azulene_job.qm_predicates,
        parameter_classes=azulene_job.parameter_classes,
    )
    fit_cases = (
        SemiexperimentalBenchmarkCase(
            "water",
            SemiexperimentalFitRequest(
                repo / "examples/semiexp/water/parent.xyz",
                read_observations(repo / "examples/semiexp/water/isotopologues.toml"),
            ),
        ),
        SemiexperimentalBenchmarkCase(
            "azulene_msr_subset",
            azulene_request,
            max_iter=azulene_job.max_iter,
            step=azulene_job.step,
            damping=azulene_job.damping,
            max_step=azulene_job.max_step,
            prune_condition=azulene_job.prune_condition,
        ),
    )
    fit_rows = []
    fit_errors = []
    for case in fit_cases:
        try:
            fit_rows.extend(run_semiexperimental_benchmark((case,), outdir=outdir / "fits", max_iter=args.max_iter))
        except Exception as exc:
            fit_errors.append((case.label, str(exc).replace("\n", " ")))
    (outdir / "semiexp_fit_benchmarks.csv").write_text(benchmark_csv(fit_rows), encoding="utf-8")
    with (outdir / "semiexp_fit_benchmark_errors.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "message"])
        writer.writerows(fit_errors)
    print(outdir / "semiexp_gic_benchmarks.csv")
    print(outdir / "semiexp_fit_benchmarks.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
