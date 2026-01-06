# paper/plot_figure3_confusion_matrices.py
"""
Generate paper Figure 3: aggregated confusion matrices for best configuration per task.

Reads:
  - performance_summary.csv
  - confusion_matrices_aggregated.csv

Writes (to --out_dir/figures):
  - Figure3_confusion_matrices.pdf
  - Figure3_confusion_matrices.jpg  (dpi=300)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


TASK_ORDER = ["CN_vs_AD", "CN_vs_EMCI", "EMCI_vs_AD"]


def _nice_task_label(task: str) -> str:
    return task.replace("_vs_", " vs ")


def _pick_best_configs(df_summary: pd.DataFrame) -> pd.DataFrame:
    """
    Pick best (domain, metric) per task based on max mean AUC.
    """
    d = df_summary[df_summary["measure"] == "auc"].copy()
    d["domain"] = d["domain"].astype(str).str.lower()
    d["metric"] = d["metric"].astype(str).str.lower()
    d = d[d["task"].isin(TASK_ORDER)]

    best_rows = []
    for t in TASK_ORDER:
        g = d[d["task"] == t].copy()
        g = g.sort_values(["mean"], ascending=False)
        best_rows.append(g.iloc[0][["task", "domain", "metric", "mean"]].to_dict())
    return pd.DataFrame(best_rows)


def run(args: argparse.Namespace) -> None:
    results_dir = Path(args.results_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    sum_path = results_dir / "performance_summary.csv"
    cm_path = results_dir / "confusion_matrices_aggregated.csv"
    if not sum_path.exists() or not cm_path.exists():
        raise FileNotFoundError("Missing result CSVs. Run scripts/nested_cv_train_eval.py first.")

    df_sum = pd.read_csv(sum_path)
    df_cm = pd.read_csv(cm_path)

    best = _pick_best_configs(df_sum)

    # Build figure
    fig, axes = plt.subplots(1, 3, figsize=(args.width, args.height))
    if len(axes.shape) == 0:
        axes = [axes]

    for ax, task in zip(axes, TASK_ORDER):
        row = best[best["task"] == task].iloc[0]
        domain = row["domain"]
        metric = row["metric"]
        mean_auc = float(row["mean"])

        cmr = df_cm[(df_cm["task"] == task) & (df_cm["domain"].str.lower() == domain) & (df_cm["metric"].str.lower() == metric)]
        if cmr.empty:
            raise ValueError(f"No confusion matrix found for {task} {domain}/{metric}")

        cmr = cmr.iloc[0]
        tp, fp, tn, fn = int(cmr["tp"]), int(cmr["fp"]), int(cmr["tn"]), int(cmr["fn"])
        total = int(cmr["total"]) if "total" in cmr else (tp + fp + tn + fn)

        # Matrix layout:
        # rows = True [neg, pos], cols = Pred [neg, pos]
        mat = np.array([[tn, fp],
                        [fn, tp]], dtype=int)

        im = ax.imshow(mat, aspect="equal")
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])

        neg = str(cmr["negative_class"])
        pos = str(cmr["positive_class"])
        ax.set_xticklabels([neg, pos], fontsize=11)
        ax.set_yticklabels([neg, pos], fontsize=11)

        ax.set_xlabel("Predicted", fontsize=11)
        ax.set_ylabel("True", fontsize=11)

        # Annotate cells with count + percent of total
        for (i, j), v in np.ndenumerate(mat):
            pct = 100.0 * v / total if total > 0 else 0.0
            ax.text(j, i, f"{v}\n({pct:.1f}%)", ha="center", va="center", fontsize=11)

        title = f"{_nice_task_label(task)}\n({domain.capitalize()} + {metric.capitalize()}, mean AUC={mean_auc:.3f})"
        if args.hide_titles:
            title = _nice_task_label(task)
        ax.set_title(title, fontsize=12)

    plt.tight_layout()

    pdf_path = fig_dir / "Figure3_confusion_matrices.pdf"
    jpg_path = fig_dir / "Figure3_confusion_matrices.jpg"
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.savefig(jpg_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Wrote: {pdf_path}")
    print(f"[OK] Wrote: {jpg_path}")


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot Figure 3 (aggregated confusion matrices).")
    p.add_argument("--results_dir", type=str, default="outputs/paper_results")
    p.add_argument("--out_dir", type=str, default="outputs/paper_outputs")
    p.add_argument("--width", type=float, default=14.5)
    p.add_argument("--height", type=float, default=4.6)
    p.add_argument("--hide_titles", action="store_true", help="Hide config/AUC in subplot titles.")
    return p


if __name__ == "__main__":
    run(build_argparser().parse_args())
