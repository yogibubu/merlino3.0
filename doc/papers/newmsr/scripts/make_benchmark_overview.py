#!/usr/bin/env python3
"""Build the model-generation audit figure for the manuscript."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "generated" / "paper_benchmark_summary.csv"
OUT_PATH = ROOT / "figures" / "benchmark_overview.pdf"

LABELS = {
    "glycolaldehyde": "Glycol.",
    "cyclopentadiene": "Cyclopent.",
    "nitrobenzene": "Nitrobenz.",
    "p-EBN": "p-EBN",
    "azulene": "Azulene",
    "norcamphor": "Norcamphor",
    "glycine_I": "Gly-I",
    "glycine_II": "Gly-II",
}

COLORS = {
    "primitive": "#4A5561",
    "active": "#245C73",
    "effective": "#2F8F83",
    "constraints": "#C88A2D",
    "ink": "#263238",
    "grid": "#D8DEE3",
    "pale": "#F4F7F8",
}


def read_rows() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    rows = read_rows()
    labels = [LABELS.get(row["system"], row["system"].replace("_", " ")) for row in rows]
    x = list(range(len(labels)))

    final_gics = [int(row["final_gics"]) for row in rows]
    active_gics = [int(row["totally_symmetric_gics"]) for row in rows]
    variables = [int(row["effective_parameters"]) for row in rows]
    constraints = [int(row["primitive_constraints"]) for row in rows]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "legend.fontsize": 7.3,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, (ax_model, ax_constraints) = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.35),
        gridspec_kw={"width_ratios": [2.25, 1.0]},
    )
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.28, top=0.60, wspace=0.23)

    box = FancyBboxPatch(
        (0.045, 0.755),
        0.91,
        0.205,
        transform=fig.transFigure,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=0.7,
        edgecolor=COLORS["grid"],
        facecolor=COLORS["pale"],
        zorder=0,
    )
    fig.add_artist(box)
    fig.text(
        0.17,
        0.902,
        "Manual model\nformulation",
        color=COLORS["ink"],
        fontsize=9.0,
        weight="bold",
        ha="center",
        va="center",
        linespacing=0.95,
    )
    fig.text(
        0.50,
        0.838,
        "topology + symmetry + coverage",
        color=COLORS["primitive"],
        fontsize=8.5,
        ha="center",
        va="center",
    )
    fig.text(
        0.81,
        0.902,
        "algorithmic model",
        color=COLORS["ink"],
        fontsize=9.4,
        weight="bold",
        ha="center",
        va="center",
    )
    fig.text(
        0.17,
        0.789,
        "0 molecule-specific\ncoordinate edits",
        color=COLORS["constraints"],
        fontsize=7.8,
        weight="bold",
        ha="center",
        va="center",
        linespacing=0.95,
    )
    fig.text(
        0.81,
        0.789,
        "primitive -> active ->\neffective variables",
        color=COLORS["effective"],
        fontsize=7.8,
        ha="center",
        va="center",
        linespacing=0.95,
    )
    fig.add_artist(
        FancyArrowPatch(
            (0.30, 0.888),
            (0.69, 0.888),
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.15,
            color=COLORS["primitive"],
        )
    )

    width = 0.24
    ax_model.bar([v - width for v in x], final_gics, width, color=COLORS["primitive"], label="Generated GICs")
    ax_model.bar(x, active_gics, width, color=COLORS["active"], label="Symmetry active")
    ax_model.bar([v + width for v in x], variables, width, color=COLORS["effective"], label="Fit variables")
    ax_model.set_title("Generated model size")
    ax_model.set_ylabel("Count")
    ax_model.set_xticks(x)
    ax_model.set_xticklabels(labels, rotation=35, ha="right")
    ax_model.grid(axis="y", color=COLORS["grid"], linewidth=0.6)
    ax_model.set_axisbelow(True)
    handles, labels_for_legend = ax_model.get_legend_handles_labels()
    fig.legend(
        handles,
        labels_for_legend,
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.37, 0.715),
        columnspacing=1.15,
        handlelength=1.25,
        handletextpad=0.45,
    )

    ax_constraints.bar(x, constraints, width=0.56, color=COLORS["constraints"])
    ax_constraints.set_title("Generated constraints")
    ax_constraints.set_xticks(x)
    ax_constraints.set_xticklabels(labels, rotation=35, ha="right")
    ax_constraints.grid(axis="y", color=COLORS["grid"], linewidth=0.6)
    ax_constraints.set_axisbelow(True)
    ax_constraints.set_ylim(0, max(constraints) * 1.22 if max(constraints) else 1)
    for i, value in enumerate(constraints):
        if value:
            ax_constraints.text(i, value + max(constraints) * 0.035, str(value), ha="center", va="bottom", fontsize=7)

    for ax in (ax_model, ax_constraints):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(COLORS["grid"])
        ax.spines["bottom"].set_color(COLORS["grid"])
        ax.tick_params(axis="both", colors=COLORS["ink"], width=0.6)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, bbox_inches="tight")


if __name__ == "__main__":
    main()
