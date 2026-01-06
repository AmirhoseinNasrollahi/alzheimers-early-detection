# scripts/nested_cv_train_eval.py
"""
Repeated stratified nested CV evaluation on similarity-score features.

Inputs:
  - similarity_features.csv with columns:
      subject_id, diagnosis, roi, domain, metric, similarity_score

Outputs (written to --out_dir):
  - performance_outerfolds.csv
  - performance_summary.csv
  - confusion_matrices_aggregated.csv

This script reproduces paper-ready evaluation artifacts WITHOUT requiring raw MRI.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score


@dataclass(frozen=True)
class BinaryTask:
    name: str
    negative_class: str
    positive_class: str


TASKS: List[BinaryTask] = [
    BinaryTask("CN_vs_AD", negative_class="CN", positive_class="AD"),
    BinaryTask("CN_vs_EMCI", negative_class="CN", positive_class="EMCI"),
    BinaryTask("EMCI_vs_AD", negative_class="EMCI", positive_class="AD"),
]


def _balanced_accuracy(tp: int, fp: int, tn: int, fn: int) -> float:
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    return 0.5 * (sens + spec)


def _precision(tp: int, fp: int) -> float:
    return tp / (tp + fp) if (tp + fp) > 0 else np.nan


def _f1(tp: int, fp: int, fn: int) -> float:
    p = _precision(tp, fp)
    r = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    if np.isnan(p) or np.isnan(r) or (p + r) == 0:
        return np.nan
    return 2 * p * r / (p + r)


def _confusion_counts(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[int, int, int, int]:
    # y_true/y_pred in {0,1}
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    return tp, fp, tn, fn


def _select_threshold_inner_cv(
    scores: np.ndarray,
    y: np.ndarray,
    inner_splits: int,
    seed: int,
    n_thresholds: int = 201,
) -> float:
    """
    Select a threshold on the ORIGINAL similarity score (scores),
    with rule: predict positive (disease) if score <= threshold.

    Objective: maximize mean balanced accuracy across inner folds.
    """
    if len(np.unique(scores)) < 2:
        # Degenerate case
        return float(np.median(scores))

    # Candidate thresholds: quantiles across training scores
    qs = np.linspace(0.0, 1.0, n_thresholds)
    candidates = np.quantile(scores, qs)
    candidates = np.unique(candidates)

    inner_cv = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=seed)

    best_thr = float(candidates[len(candidates) // 2])
    best_score = -np.inf

    for thr in candidates:
        fold_scores = []
        for tr_idx, va_idx in inner_cv.split(scores, y):
            s_tr, s_va = scores[tr_idx], scores[va_idx]
            y_va = y[va_idx]

            y_pred_va = (s_va <= thr).astype(int)
            tp, fp, tn, fn = _confusion_counts(y_va, y_pred_va)
            fold_scores.append(_balanced_accuracy(tp, fp, tn, fn))

        mean_bacc = float(np.nanmean(fold_scores))
        if mean_bacc > best_score:
            best_score = mean_bacc
            best_thr = float(thr)

    return best_thr


def _mean_ci95(values: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Mean, SD, and 95% CI via t-interval across outer evaluations.
    """
    values = np.asarray(values, dtype=float)
    n = int(np.sum(~np.isnan(values)))
    mean = float(np.nanmean(values))
    sd = float(np.nanstd(values, ddof=1)) if n > 1 else float("nan")
    if n <= 1:
        return mean, sd, float("nan"), float("nan")

    # t critical
    try:
        from scipy.stats import t
        tcrit = float(t.ppf(0.975, df=n - 1))
    except Exception:
        # fallback approx
        tcrit = 1.96

    half = tcrit * sd / np.sqrt(n)
    return mean, sd, mean - half, mean + half


def run(args: argparse.Namespace) -> None:
    sim_path = Path(args.similarity_csv).expanduser().resolve()
    if not sim_path.exists():
        raise FileNotFoundError(f"similarity CSV not found: {sim_path}")

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(sim_path)

    required_cols = {"subject_id", "diagnosis", "roi", "domain", "metric", "similarity_score"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {sim_path.name}: {sorted(missing)}")

    # Filter ROI if provided
    if args.roi:
        df = df[df["roi"] == args.roi].copy()
        if df.empty:
            raise ValueError(f"No rows found for roi='{args.roi}' in {sim_path}")

    # Normalize labels (defensive)
    df["diagnosis"] = df["diagnosis"].astype(str).str.strip()
    df["domain"] = df["domain"].astype(str).str.strip().str.lower()
    df["metric"] = df["metric"].astype(str).str.strip().str.lower()

    configs = df[["domain", "metric"]].drop_duplicates().sort_values(["domain", "metric"]).to_records(index=False)
    configs = [(str(d), str(m)) for d, m in configs]

    outer_cv = RepeatedStratifiedKFold(
        n_splits=args.outer_splits,
        n_repeats=args.outer_repeats,
        random_state=args.seed,
    )

    outer_rows: List[Dict] = []
    conf_rows: List[Dict] = []

    for task in TASKS:
        # Task subset
        df_task = df[df["diagnosis"].isin([task.negative_class, task.positive_class])].copy()
        if df_task.empty:
            raise ValueError(f"No data for task {task.name} with classes {task.negative_class}/{task.positive_class}")

        # y: 1 = positive (disease), 0 = negative
        df_task["y"] = (df_task["diagnosis"] == task.positive_class).astype(int)

        for domain, metric in configs:
            df_cfg = df_task[(df_task["domain"] == domain) & (df_task["metric"] == metric)].copy()
            if df_cfg.empty:
                continue

            # IMPORTANT: classification uses ORIGINAL similarity_score (higher for CN, lower for disease)
            scores = df_cfg["similarity_score"].to_numpy(dtype=float)
            y = df_cfg["y"].to_numpy(dtype=int)

            # AUC uses inverted score so that higher = more likely disease
            inv_scores = -scores

            split_idx = 0
            for tr_idx, te_idx in outer_cv.split(scores, y):
                split_idx += 1
                repeat = (split_idx - 1) // args.outer_splits + 1
                outer_fold = (split_idx - 1) % args.outer_splits + 1

                s_tr, s_te = scores[tr_idx], scores[te_idx]
                y_tr, y_te = y[tr_idx], y[te_idx]
                inv_te = inv_scores[te_idx]

                # Inner loop: select threshold on s_tr
                thr = _select_threshold_inner_cv(
                    scores=s_tr,
                    y=y_tr,
                    inner_splits=args.inner_splits,
                    seed=args.seed + repeat * 100 + outer_fold,
                    n_thresholds=args.n_thresholds,
                )

                # Predict: disease if similarity <= threshold
                y_pred = (s_te <= thr).astype(int)

                tp, fp, tn, fn = _confusion_counts(y_te, y_pred)

                auc = float(roc_auc_score(y_te, inv_te)) if len(np.unique(y_te)) == 2 else float("nan")
                sens = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
                spec = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
                bacc = _balanced_accuracy(tp, fp, tn, fn)

                outer_rows.append(
                    dict(
                        task=task.name,
                        domain=domain,
                        metric=metric,
                        repeat=repeat,
                        outer_fold=outer_fold,
                        selected_threshold=float(thr),
                        confidence_level=float(args.confidence_level),
                        auc=auc,
                        balanced_accuracy=float(bacc),
                        sensitivity=float(sens),
                        specificity=float(spec),
                    )
                )

                conf_rows.append(
                    dict(
                        task=task.name,
                        domain=domain,
                        metric=metric,
                        positive_class=task.positive_class,
                        negative_class=task.negative_class,
                        tp=tp,
                        fp=fp,
                        tn=tn,
                        fn=fn,
                    )
                )

    df_outer = pd.DataFrame(outer_rows)
    df_conf = pd.DataFrame(conf_rows)

    # Aggregate confusion matrices across outer folds & repetitions
    df_conf_agg = (
        df_conf.groupby(["task", "domain", "metric", "positive_class", "negative_class"], as_index=False)[["tp", "fp", "tn", "fn"]]
        .sum()
    )
    df_conf_agg["total"] = df_conf_agg["tp"] + df_conf_agg["fp"] + df_conf_agg["tn"] + df_conf_agg["fn"]

    # Summary stats (mean/sd/ci95) across outer evaluations
    summary_rows: List[Dict] = []
    for (task, domain, metric), g in df_outer.groupby(["task", "domain", "metric"]):
        n = int(len(g))
        for measure in ["auc", "balanced_accuracy", "sensitivity", "specificity"]:
            mean, sd, lo, hi = _mean_ci95(g[measure].to_numpy(dtype=float))
            summary_rows.append(
                dict(
                    task=task,
                    domain=domain,
                    metric=metric,
                    measure=measure,
                    mean=mean,
                    sd=sd,
                    ci95_low=lo,
                    ci95_high=hi,
                    n_outer_evals=n,
                )
            )
    df_summary = pd.DataFrame(summary_rows)

    # Write outputs
    out_outer = out_dir / "performance_outerfolds.csv"
    out_summary = out_dir / "performance_summary.csv"
    out_conf = out_dir / "confusion_matrices_aggregated.csv"

    df_outer.to_csv(out_outer, index=False)
    df_summary.to_csv(out_summary, index=False)
    df_conf_agg.to_csv(out_conf, index=False)

    print(f"[OK] Wrote: {out_outer}")
    print(f"[OK] Wrote: {out_summary}")
    print(f"[OK] Wrote: {out_conf}")


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Repeated nested CV evaluation on similarity-score features.")
    p.add_argument(
        "--similarity_csv",
        type=str,
        default="data/derived/similarity_features.csv",
        help="Path to similarity_features.csv (derived features).",
    )
    p.add_argument(
        "--roi",
        type=str,
        default="hippocampus_bilateral",
        help="ROI name to filter (set empty string to disable filtering).",
    )
    p.add_argument(
        "--out_dir",
        type=str,
        default="outputs/paper_results",
        help="Output directory for result CSV files.",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--outer_splits", type=int, default=5)
    p.add_argument("--outer_repeats", type=int, default=5)
    p.add_argument("--inner_splits", type=int, default=5)
    p.add_argument("--n_thresholds", type=int, default=201, help="Number of quantile-based threshold candidates.")
    p.add_argument("--confidence_level", type=float, default=0.95, help="Stored in CSV for reporting.")
    return p


if __name__ == "__main__":
    args = build_argparser().parse_args()
    # Allow disabling ROI filtering by passing --roi ""
    if args.roi is not None and str(args.roi).strip() == "":
        args.roi = None
    run(args)
