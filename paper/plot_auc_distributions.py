# paper/plot_figure_auc_distributions.py
"""
Generate paper Figure : AUC distributions across repeated nested CV.

Reads:
  - performance_outerfolds.csv

Writes (to --out_dir/figures):
  - Figure2_auc_distributions.pdf
  - Figure2_auc_distributions.jpg  (dpi=300)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


TASK_ORDER = ["CN_vs_AD", "CN_vs_EMCI", "EMCI_vs_AD"]
METRIC_ORDER = ["pearson", "spearman", "kendall", "distance"]
DOMAIN_ORDER = ["spatial", "frequency"]


def _nice_task_label(task: str) -> str:
    return task.replace("_vs_", " vs ").replace("CN", "CN").replace("EMCI", "EMCI").replace("AD", "AD")


def run(args: argparse.Namespace) -> None:
    results_dir = Path(args.results_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    csv_path = results_dir / "performance_outerfolds.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing {csv_path}. Run scripts/nested_cv_train_eval.py first.")

    df = pd.read_csv(csv_path)
    df["task"] = df["task"].astype(str)
    df["domain"] = df["domain"].astype(str).str.lower()
    df["metric"] = df["metric"].astype(str).str.lower()

    # Keep only AUC and known tasks
    df = df[df["task"].isin(TASK_ORDER)].copy()

    # Figure sizing: bigger = more journal-friendly
    fig_w = args.width
    fig_h = args.height
    plt.figure(figsize=(fig_w, fig_h))

    # Build x positions: blocks per task, metrics within
    block_gap = 1.8
    metric_gap = 1.0
    domain_offset = 0.18  # side-by-side within metric

    x_positions = {}
    xticks = []
    xticklabels = []

    base = 0.0
    for t in TASK_ORDER:
        for mi, m in enumerate(METRIC_ORDER):
            x = base + mi * metric_gap
            x_positions[(t, m)] = x
            xticks.append(x)
            xticklabels.append(m.capitalize() if m != "distance" else "Distance corr.")
        base += (len(METRIC_ORDER) - 1) * metric_gap + block_gap

    ax = plt.gca()

    # Violin plots
    for t in TASK_ORDER:
        for m in METRIC_ORDER:
            x_center = x_positions[(t, m)]
            for di, d in enumerate(DOMAIN_ORDER):
                vals = df[(df["task"] == t) & (df["metric"] == m) & (df["domain"] == d)]["auc"].to_numpy()
                if vals.size == 0:
                    continue

                x = x_center + (-domain_offset if d == "spatial" else domain_offset)

                parts = ax.violinplot(
                    dataset=[vals],
                    positions=[x],
                    widths=0.28,
                    showmeans=False,
                    showmedians=False,
                    showextrema=False,
                )

                # Make violins clean (no custom colors; journal-friendly defaults)
                for pc in parts["bodies"]:
                    pc.set_alpha(0.20)

                # Box overlay
                bp = ax.boxplot(
                    [vals],
                    positions=[x],
                    widths=0.10,
                    patch_artist=True,
                    showfliers=False,
                    medianprops=dict(linewidth=1.5),
                    boxprops=dict(linewidth=1.0),
                    whiskerprops=dict(linewidth=1.0),
                    capprops=dict(linewidth=1.0),
                )
                for box in bp["boxes"]:
                    box.set_alpha(0.25)

                # Jittered points
                rng = np.random.default_rng(123)
                jitter = rng.normal(0, 0.03, size=vals.size)
                marker = "o" if d == "spatial" else "^"
                ax.scatter(np.full(vals.size, x) + jitter, vals, s=18, alpha=0.85, marker=marker)

    # Task labels under x-axis
    base = 0.0
    for t in TASK_ORDER:
        mid = base + (len(METRIC_ORDER) - 1) * metric_gap / 2
        ax.text(mid, -0.12, _nice_task_label(t), transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=12)
        base += (len(METRIC_ORDER) - 1) * metric_gap + block_gap

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=11)
    ax.set_ylabel("AUC (outer-fold evaluations)", fontsize=12)

    ax.set_ylim(args.ymin, args.ymax)
    ax.grid(True, axis="y", alpha=0.25)

    # Legend (marker-based)
    ax.scatter([], [], marker="o", label="Spatial")
    ax.scatter([], [], marker="^", label="Frequency")
    ax.legend(frameon=True, loc="upper left")

    plt.tight_layout()

    pdf_path = fig_dir / "Figure2_auc_distributions.pdf"
    jpg_path = fig_dir / "Figure2_auc_distributions.jpg"
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.savefig(jpg_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Wrote: {pdf_path}")
    print(f"[OK] Wrote: {jpg_path}")


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot Figure 2 (AUC distributions).")
    p.add_argument("--results_dir", type=str, default="outputs/paper_results")
    p.add_argument("--out_dir", type=str, default="outputs/paper_outputs")
    p.add_argument("--width", type=float, default=16.0)
    p.add_argument("--height", type=float, default=4.8)
    p.add_argument("--ymin", type=float, default=0.50)
    p.add_argument("--ymax", type=float, default=1.00)
    return p


if __name__ == "__main__":
    run(build_argparser().parse_args())
