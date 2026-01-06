#!/usr/bin/env bash
set -euo pipefail

# paper/reproduce_all.sh
# Reproduce paper-ready results, tables, and figures from derived similarity features.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SIM_CSV="${ROOT_DIR}/data/derived/similarity_features.csv"
RESULTS_DIR="${ROOT_DIR}/outputs/paper_results"
OUT_DIR="${ROOT_DIR}/outputs/paper_outputs"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${OUT_DIR}"

echo "[1/4] Running nested CV evaluation..."
python "${ROOT_DIR}/scripts/nested_cv_train_eval.py" \
  --similarity_csv "${SIM_CSV}" \
  --out_dir "${RESULTS_DIR}" \
  --roi "hippocampus_bilateral" \
  --outer_splits 5 --outer_repeats 5 --inner_splits 5 \
  --seed 42

echo "[2/4] Building tables (Table2 & Table3)..."
python "${ROOT_DIR}/paper/build_tables.py" \
  --results_dir "${RESULTS_DIR}" \
  --out_dir "${OUT_DIR}"

echo "[3/4] Plotting Figure2 (AUC distributions)..."
python "${ROOT_DIR}/paper/plot_figure2_auc_distributions.py" \
  --results_dir "${RESULTS_DIR}" \
  --out_dir "${OUT_DIR}"

echo "[4/4] Plotting Figure3 (confusion matrices)..."
python "${ROOT_DIR}/paper/plot_figure3_confusion_matrices.py" \
  --results_dir "${RESULTS_DIR}" \
  --out_dir "${OUT_DIR}"

echo "[DONE] Outputs written to:"
echo "  Results: ${RESULTS_DIR}"
echo "  Paper outputs (figures/tables): ${OUT_DIR}"
