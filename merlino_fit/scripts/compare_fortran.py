"""Placeholder for regression vs Fortran reference.

Expected usage:
    python scripts/compare_fortran.py --data data.npz --fortran-out ref_terms.txt

This script will compare fitted terms vs a reference file from the legacy
Fortran pipeline (mkprim/mksalc-based). The Fortran code is used only as a
numerical reference.
"""

import argparse

import numpy as np


def load_terms(path):
    terms = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            idxs = tuple(int(p) for p in parts[:-1])
            val = float(parts[-1])
            terms[idxs] = val
    return terms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-out", required=True)
    ap.add_argument("--fortran-out", required=True)
    args = ap.parse_args()

    t_fit = load_terms(args.fit_out)
    t_ref = load_terms(args.fortran_out)

    common = set(t_fit.keys()) & set(t_ref.keys())
    if not common:
        raise SystemExit("No common terms between fit and Fortran outputs")

    diffs = [abs(t_fit[k] - t_ref[k]) for k in common]
    print(f"Common terms: {len(common)}")
    print(f"Max diff: {max(diffs):.6e}")
    print(f"RMS diff: {(np.mean(np.square(diffs))**0.5):.6e}")


if __name__ == "__main__":
    main()
