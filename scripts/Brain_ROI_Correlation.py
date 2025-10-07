#!/usr/bin/env python3
"""
Brain ROI Correlation & FFT CLI
================================

This tool computes **voxelwise correlations** between 3D case ROIs and atlas ROIs
in the spatial domain and in the frequency domain (via 3D FFT).

For each case:
- Loads the case ROI volume and the corresponding atlas ROI volume
- Computes Pearson, Spearman, Kendall, and Distance Correlation (via `dcor` if available; otherwise a fallback) on **flattened 3D arrays**
- Saves results to two CSVs in a per‑case folder:
  - `<CASE_NAME>+correlation.csv` (spatial domain)
  - `<CASE_NAME>+FFT-correlation.csv` (frequency domain on |FFTN| magnitudes)
- Saves the frequency‑domain magnitude image as `<CASE_NAME>+FFTN.nii.gz`

Designed for **3D** medical images (e.g., NIfTI `.nii/.nii.gz`, NRRD `.nrrd`, MHA/MHD).

Author: You 🚀
License: MIT (customize if needed)
"""
from __future__ import annotations

import argparse
import os
import sys
import math
import glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# --- Optional imports ---
try:
    import dcor  # true distance correlation
    _HAS_DCOR = True
except Exception:
    _HAS_DCOR = False

try:
    from tqdm import tqdm
    _HAS_TQDM = True
except Exception:
    _HAS_TQDM = False

import SimpleITK as sitk
from scipy.stats import pearsonr, spearmanr, kendalltau


# --------------------------------------
# Utilities
# --------------------------------------

SUPPORTED_EXTS = {".nii", ".nii.gz", ".nrrd", ".mha", ".mhd"}


def is_image_file(p: Path) -> bool:
    """Return True if path has a supported medical image extension."""
    name = p.name.lower()
    return any(name.endswith(ext) for ext in SUPPORTED_EXTS)


def read_image(path: Path) -> Tuple[np.ndarray, sitk.Image]:
    """Read medical image via SimpleITK -> (numpy array [z,y,x], original sitk.Image)."""
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img)
    # Ensure float for downstream computations
    if not np.issubdtype(arr.dtype, np.floating):
        arr = arr.astype(np.float32, copy=False)
    return arr, img


def write_like(reference: sitk.Image, array: np.ndarray, out_path: Path) -> None:
    """Write numpy array back as image using reference metadata (spacing, origin, direction)."""
    sitk_img = sitk.GetImageFromArray(array)
    sitk_img.CopyInformation(reference)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(sitk_img, str(out_path))


def zscore(x: np.ndarray) -> np.ndarray:
    """Z-score normalize array (over all voxels)."""
    m = np.mean(x)
    s = np.std(x)
    if s == 0:
        return np.zeros_like(x, dtype=np.float32)
    return (x - m) / s


def distance_correlation_fallback(x: np.ndarray, y: np.ndarray) -> float:
    """
    Lightweight distance correlation approximation (no full pairwise distance matrices).
    NOTE: Prefer dcor.distance_correlation for correctness when available.
    """
    # mean-centered absolute distances to the mean (not full pairwise distances)
    ax = np.abs(x - x.mean())
    ay = np.abs(y - y.mean())
    dcov = np.mean(ax * ay)
    dvarx = np.mean(ax**2)
    dvary = np.mean(ay**2)
    denom = math.sqrt(dvarx * dvary) if dvarx > 0 and dvary > 0 else 0.0
    return float(dcov / denom) if denom != 0 else 0.0


def distance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Distance correlation via dcor if available; else fallback approximation."""
    if _HAS_DCOR:
        try:
            return float(dcor.distance_correlation(x, y))
        except Exception:
            pass
    return distance_correlation_fallback(x, y)


def compute_correlations(a: np.ndarray, b: np.ndarray) -> Dict[str, float]:
    """Compute Pearson, Spearman, Kendall, Distance on flattened arrays."""
    af = a.reshape(-1)
    bf = b.reshape(-1)

    # Guard: drop NaNs/inf if present
    mask = np.isfinite(af) & np.isfinite(bf)
    af = af[mask]
    bf = bf[mask]

    # Stats
    pr, _ = pearsonr(af, bf)
    sr, _ = spearmanr(af, bf)
    kt, _ = kendalltau(af, bf)
    dc = distance_correlation(af, bf)

    return {"Pearson": float(pr), "Spearman": float(sr), "Kendall": float(kt), "Distance": float(dc)}


def fft_magnitude(vol: np.ndarray, shift: bool = False) -> np.ndarray:
    """Compute |FFTN| magnitude of a 3D volume. Optionally fftshift for nicer anatomy of spectrum."""
    F = np.fft.fftn(vol, axes=(0, 1, 2))
    if shift:
        F = np.fft.fftshift(F, axes=(0, 1, 2))
    mag = np.abs(F)
    return mag.astype(np.float32, copy=False)


def find_images(path: Path) -> List[Path]:
    """Find supported image files under path (single file or directory)."""
    if path.is_file():
        return [path] if is_image_file(path) else []
    files = []
    for ext in [".nii", ".nii.gz", ".nrrd", ".mha", ".mhd"]:
        files.extend(path.rglob(f"*{ext}"))
    return sorted(files)


def stem_key(p: Path) -> str:
    """Robust stem: remove double extensions like .nii.gz and keep core name."""
    name = p.name
    for ext in (".nii.gz", ".nii", ".nrrd", ".mha", ".mhd"):
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return p.stem


def build_atlas_map(atlas_files: List[Path]) -> Dict[str, Path]:
    """
    Build a mapping from stem -> atlas file.
    If only ONE atlas file exists, map it under the special key '*ALL*' to reuse for every case.
    """
    if len(atlas_files) == 1:
        return {"*ALL*": atlas_files[0]}
    d: Dict[str, Path] = {}
    for f in atlas_files:
        d[stem_key(f).lower()] = f
    return d


def match_atlas_for_case(case_path: Path, atlas_index: Dict[str, Path]) -> Optional[Path]:
    """Find the best atlas file for the given case by stem matching; fallback to single atlas."""
    if "*ALL*" in atlas_index:
        return atlas_index["*ALL*"]
    cstem = stem_key(case_path).lower()
    # exact stem match first
    if cstem in atlas_index:
        return atlas_index[cstem]
    # fuzzy: substring match
    for k, p in atlas_index.items():
        if k in cstem or cstem in k:
            return p
    return None


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def save_csv(d: Dict[str, float], out_csv: Path) -> None:
    df = pd.DataFrame({"Metric": list(d.keys()), "Correlation": list(d.values())})
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)


# --------------------------------------
# Main pipeline
# --------------------------------------

def process_case(case_path: Path,
                 atlas_path: Path,
                 out_root: Path,
                 normalize: bool = False,
                 mask_nonzero_union: bool = False,
                 fft_shift: bool = False) -> None:
    """
    Process a single case file against a chosen atlas file.
    Saves:
      - <name>/<name>+correlation.csv
      - <name>/<name>+FFTN.nii.gz
      - <name>/<name>+FFT-correlation.csv
    """
    case_name = stem_key(case_path)
    case_out_dir = out_root / case_name
    ensure_dir(case_out_dir)

    # --- Load images ---
    case_vol, case_img = read_image(case_path)
    atlas_vol, atlas_img = read_image(atlas_path)

    # Sanity: shapes must match
    if case_vol.shape != atlas_vol.shape:
        raise ValueError(f"Shape mismatch for case '{case_name}': "
                         f"case {case_vol.shape} vs atlas {atlas_vol.shape}")

    # Optional normalization
    if normalize:
        case_proc = zscore(case_vol)
        atlas_proc = zscore(atlas_vol)
    else:
        case_proc = case_vol
        atlas_proc = atlas_vol

    # Optional masking to reduce zero-background effects
    if mask_nonzero_union:
        mask = (case_proc != 0) | (atlas_proc != 0)
        case_use = case_proc[mask]
        atlas_use = atlas_proc[mask]
        # Keep full arrays for FFT
        case_full = case_proc
        atlas_full = atlas_proc
    else:
        case_use = case_proc
        atlas_use = atlas_proc
        case_full = case_proc
        atlas_full = atlas_proc

    # --- Spatial correlations ---
    spatial_corrs = compute_correlations(case_use, atlas_use)
    save_csv(spatial_corrs, case_out_dir / f"{case_name}+correlation.csv")

    # --- FFT magnitude & save ---
    case_fft_mag = fft_magnitude(case_full, shift=fft_shift)
    atlas_fft_mag = fft_magnitude(atlas_full, shift=fft_shift)
    fft_out_path = case_out_dir / f"{case_name}+FFTN.nii.gz"
    write_like(case_img, case_fft_mag, fft_out_path)

    # --- Frequency-domain correlations (on magnitudes) ---
    freq_corrs = compute_correlations(case_fft_mag, atlas_fft_mag)
    save_csv(freq_corrs, case_out_dir / f"{case_name}+FFT-correlation.csv")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Compute 3D ROI correlations (spatial & frequency) between case volumes and atlas volumes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--cases", required=True, type=Path,
                   help="Path to a **directory** containing case ROI volumes, or a single case file.")
    p.add_argument("--atlas", required=True, type=Path,
                   help="Path to a **directory** containing atlas ROI volumes, or a single atlas file.")
    p.add_argument("--out", required=True, type=Path,
                   help="Output root directory. Per-case subfolders will be created here.")
    p.add_argument("--normalize", action="store_true",
                   help="Z-score normalize volumes before correlation.")
    p.add_argument("--mask-nonzero-union", action="store_true",
                   help="Only use voxels where either (case != 0) OR (atlas != 0) for spatial correlations.")
    p.add_argument("--fft-shift", action="store_true",
                   help="Apply fftshift before taking magnitude (affects saved FFT image and frequency correlations).")
    p.add_argument("--strict-match", action="store_true",
                   help="Require exact stem name match between case and atlas when atlas has multiple files.")
    args = p.parse_args(argv)

    case_files = find_images(args.cases)
    if not case_files:
        p.error(f"No case images were found under: {args.cases}")

    atlas_files = find_images(args.atlas)
    if not atlas_files:
        p.error(f"No atlas images were found under: {args.atlas}")

    atlas_index = build_atlas_map(atlas_files)

    it = case_files
    if _HAS_TQDM:
        it = tqdm(case_files, desc="Processing cases", unit="case")

    for case_path in it:
        atlas_path = match_atlas_for_case(case_path, atlas_index)
        if atlas_path is None:
            if args.strict_match:
                raise SystemExit(f"[strict-match] Could not find matching atlas for case: {case_path.name}")
            # fallback: use the first atlas file
            atlas_path = atlas_files[0]

        try:
            process_case(case_path, atlas_path, args.out,
                         normalize=args.normalize,
                         mask_nonzero_union=args.mask_nonzero_union,
                         fft_shift=args.fft_shift)
        except Exception as e:
            print(f"[ERROR] Case '{case_path.name}': {e}", file=sys.stderr)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
