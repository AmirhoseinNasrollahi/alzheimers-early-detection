#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FFT & Spectral Analysis Pipeline for Brain MRI Volumes (NIfTI)
- CLI-friendly and GitHub-ready
- Mirrors the input folder structure in the output directory
- Two stages:
    1) fft:     compute 3D FFT power spectrum volumes from spatial-domain NIfTI
    2) analyze: perform spectral analysis on FFT-domain volumes
- Option 'both' runs the two stages sequentially.

Author: Your Name
License: MIT
"""

import os
import sys
import argparse
import json
import math
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------- Defaults -----------------------------

N_BINS_PROFILE = 256     # radial bins for radial average profile
N_BANDS = 8              # number of radial bands for band energy
N_BOOT = 2000            # bootstrap resamples
CI_LEVEL = 0.95
RANDOM_SEED = 42

VALID_EXTS = (".nii", ".nii.gz")


# ----------------------------- Utils -----------------------------

def find_niis(root: Path, pattern: str = "*.nii*") -> List[Path]:
    """Recursively list NIfTI files under root matching pattern."""
    return sorted([p for p in root.rglob(pattern) if p.suffix in (".nii",) or "".join(p.suffixes).endswith(".nii.gz")])


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def load_nii(path: Path) -> Tuple[np.ndarray, nib.Nifti1Image]:
    """Load NIfTI and return ndarray [Z,Y,X] and the original image object (for affine/header reuse)."""
    img = nib.load(str(path))
    arr = img.get_fdata().astype(np.float32)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

    # Heuristic: Many NIfTI volumes are [X,Y,Z]. Convert to [Z,Y,X] for axial slicing.
    # If last axis is the largest, we treat it as Z and transpose.
    if arr.ndim != 3:
        raise ValueError(f"Expected 3D volume at {path}, got shape {arr.shape}")
    if arr.shape[2] >= arr.shape[0] and arr.shape[2] >= arr.shape[1]:
        arr_zyx = np.transpose(arr, (2, 1, 0))  # [Z, Y, X]
    else:
        arr_zyx = arr

    return arr_zyx, img


def save_nii_like(out_path: Path, arr_zyx: np.ndarray, ref_img: nib.Nifti1Image) -> None:
    """
    Save a [Z,Y,X] array back to NIfTI using ref_img affine/header.
    We transpose back to [X,Y,Z] to keep consistency with common viewers.
    """
    ensure_parent(out_path)
    arr_xyz = np.transpose(arr_zyx, (2, 1, 0))
    new_img = nib.Nifti1Image(arr_xyz, affine=ref_img.affine, header=ref_img.header)
    nib.save(new_img, str(out_path))


def fft_power_spectrum_3d(vol_zyx: np.ndarray) -> np.ndarray:
    """
    Compute 3D FFT power spectrum (magnitude squared), shifted to center low freq.
    Returns float32 array same shape as vol_zyx.
    """
    # Real-to-complex FFT:
    F = np.fft.fftn(vol_zyx, axes=(0, 1, 2))
    F = np.fft.fftshift(F, axes=(0, 1, 2))
    P = np.abs(F) ** 2
    # Optional dynamic range compression:
    P = np.log1p(P)  # log(1+|F|^2)
    return P.astype(np.float32)


def radial_profile_2d(img2d: np.ndarray, n_bins: int = N_BINS_PROFILE) -> Tuple[np.ndarray, np.ndarray]:
    """Compute radial average profile on a centered 2D image."""
    h, w = img2d.shape
    cy, cx = (h - 1) / 2.0, (w - 1) / 2.0
    y, x = np.indices((h, w))
    r = np.sqrt((y - cy) ** 2 + (x - cx) ** 2)
    r_norm = r / np.max(r) if np.max(r) > 0 else r
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.digitize(r_norm.ravel(), bins) - 1
    idx = np.clip(idx, 0, n_bins - 1)

    vals = img2d.ravel()
    sums = np.bincount(idx, weights=vals, minlength=n_bins)
    counts = np.bincount(idx, minlength=n_bins).astype(np.float32)
    means = np.where(counts > 0, sums / counts, 0.0)
    radii_centers = 0.5 * (bins[:-1] + bins[1:])
    return radii_centers, means


def per_slice_radial_profiles(vol_zyx: np.ndarray, n_bins: int = N_BINS_PROFILE) -> Tuple[np.ndarray, np.ndarray]:
    """Compute radial profiles for each axial slice. Returns r_axis [B], profiles [Z,B]."""
    z = vol_zyx.shape[0]
    profs = np.zeros((z, n_bins), dtype=np.float32)
    r_axis = None
    for k in range(z):
        r_axis, m = radial_profile_2d(vol_zyx[k], n_bins)
        profs[k] = m
    if r_axis is None:
        r_axis, _ = radial_profile_2d(vol_zyx[min(0, z-1)], n_bins)
    return r_axis, profs


def l1_norm_rows(arr2d: np.ndarray) -> np.ndarray:
    s = np.sum(np.abs(arr2d), axis=1, keepdims=True)
    s[s == 0] = 1.0
    return arr2d / s


def band_edges_centers(n_bands: int = N_BANDS):
    edges = np.linspace(0.0, 1.0, n_bands + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    labels = [f"[{edges[i]:.2f}-{edges[i+1]:.2f})" for i in range(n_bands)]
    return edges, centers, labels


def per_slice_band_energies(vol_zyx: np.ndarray, n_bands: int = N_BANDS):
    """Compute per-slice band-limited energy and fractional energy across radial bands."""
    z, h, w = vol_zyx.shape
    cy, cx = (h - 1) / 2.0, (w - 1) / 2.0
    y, x = np.indices((h, w))
    r = np.sqrt((y - cy) ** 2 + (x - cx) ** 2)
    r_norm = r / np.max(r) if np.max(r) > 0 else r
    edges, centers, labels = band_edges_centers(n_bands)
    band_vals = np.zeros((z, n_bands), dtype=np.float64)
    for i in range(n_bands):
        mask = (r_norm >= edges[i]) & (r_norm < edges[i + 1])
        if not np.any(mask):
            continue
        band_vals[:, i] = vol_zyx[:, mask].sum(axis=1)

    # Normalize per-slice to convert to fractional band energy
    row_sums = band_vals.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    band_frac = band_vals / row_sums
    return centers, labels, band_vals, band_frac


def bootstrap_mean_ci(samples_2d: np.ndarray, n_boot: int, ci: float, rng: np.random.Generator):
    """
    samples_2d: [N, D] -> N samples (slices), D features (bins or bands)
    Returns: mean[D], low[D], high[D] for the bootstrapped mean across N
    """
    N, D = samples_2d.shape
    means = np.zeros((n_boot, D), dtype=np.float64)
    for b in range(n_boot):
        idx = rng.integers(0, N, size=N, endpoint=False)
        means[b] = samples_2d[idx].mean(axis=0)
    alpha = (1.0 - ci) / 2.0
    low = np.quantile(means, alpha, axis=0)
    high = np.quantile(means, 1.0 - alpha, axis=0)
    mean = samples_2d.mean(axis=0)
    return mean, low, high


# ----------------------------- Pipeline Steps -----------------------------

def stage_fft(input_root: Path, output_root: Path, pattern: str, manifest_rows: List[Dict]) -> None:
    """
    Compute FFT power spectrum for all NIfTI volumes found under input_root.
    Save to output_root while preserving directory structure (prefix 'fft_').
    """
    files = find_niis(input_root, pattern)
    if not files:
        print(f"[fft] No NIfTI files found under: {input_root} (pattern={pattern})")
        return

    for src in files:
        rel = src.relative_to(input_root)
        name = src.stem
        # handle .nii.gz two suffixes
        if "".join(src.suffixes).endswith(".nii.gz"):
            name = src.name[:-7]  # remove .nii.gz
        out_path = output_root / rel.parent / f"fft_{name}.nii.gz"

        try:
            vol_zyx, ref_img = load_nii(src)
            ps_zyx = fft_power_spectrum_3d(vol_zyx)
            save_nii_like(out_path, ps_zyx, ref_img)
            manifest_rows.append({
                "stage": "fft",
                "input": str(src),
                "output": str(out_path),
                "status": "ok",
                "shape_z_y_x": f"{list(ps_zyx.shape)}"
            })
            print(f"[fft] Saved: {out_path}")
        except Exception as e:
            manifest_rows.append({
                "stage": "fft",
                "input": str(src),
                "output": str(out_path),
                "status": f"error: {e}",
                "shape_z_y_x": ""
            })
            print(f"[fft][ERROR] {src}: {e}", file=sys.stderr)


def stage_analyze(input_root: Path,
                  output_root: Path,
                  pattern: str,
                  n_bins: int,
                  n_bands: int,
                  n_boot: int,
                  ci: float,
                  seed: int,
                  manifest_rows: List[Dict]) -> None:
    """
    Spectral analysis for all FFT-domain volumes found under input_root.
    For each file, creates an output folder: <output_root>/<rel_dir>/<base_no_ext>/
      - radial_profile_mean_CI.csv
      - radial_profiles_slices.csv
      - band_energies_fractional_mean_CI.csv
      - band_energies_slices.csv
      - PNG plots with 95% CI
    """
    rng = np.random.default_rng(seed)
    files = find_niis(input_root, pattern)
    if not files:
        print(f"[analyze] No FFT NIfTI files found under: {input_root} (pattern={pattern})")
        return

    for src in files:
        rel = src.relative_to(input_root)
        base = src.stem
        if "".join(src.suffixes).endswith(".nii.gz"):
            base = src.name[:-7]
        # Each file gets its own output folder to avoid mixing data
        out_dir = output_root / rel.parent / base
        out_dir.mkdir(parents=True, exist_ok=True)

        try:
            vol_zyx, _ = load_nii(src)

            # Radial profiles per slice
            r_axis, prof_slices = per_slice_radial_profiles(vol_zyx, n_bins)
            prof_slices_n = l1_norm_rows(prof_slices)
            prof_mean, prof_low, prof_high = bootstrap_mean_ci(prof_slices_n, n_boot, ci, rng)

            # Band energies per slice
            band_centers, band_labels = band_edges_centers(n_bands)[1:]
            _, _, band_raw, band_frac = per_slice_band_energies(vol_zyx, n_bands)
            band_mean, band_low, band_high = bootstrap_mean_ci(band_frac, n_boot, ci, rng)

            # CSV outputs
            df_prof = pd.DataFrame({
                "radius_norm": r_axis,
                "mean": prof_mean, "low": prof_low, "high": prof_high,
            })
            df_prof.to_csv(out_dir / "radial_profile_mean_CI.csv", index=False)

            pd.DataFrame(prof_slices_n).to_csv(out_dir / "radial_profiles_slices.csv", index=False)

            df_bands = pd.DataFrame({
                "band_label": band_labels,
                "band_center": band_centers,
                "mean": band_mean, "low": band_low, "high": band_high,
            })
            df_bands.to_csv(out_dir / "band_energies_fractional_mean_CI.csv", index=False)

            pd.DataFrame(band_frac).to_csv(out_dir / "band_energies_slices.csv", index=False)

            # Plots
            plt.figure()
            plt.plot(r_axis, prof_mean, label="mean")
            plt.fill_between(r_axis, prof_low, prof_high, alpha=0.2, label=f"{int(ci*100)}% CI")
            plt.xlabel("Normalized radial spatial frequency")
            plt.ylabel("Normalized power (area = 1)")
            plt.title("Radial Average Power Spectrum (mean ± CI)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(out_dir / "radial_profiles_mean_CI.png", dpi=300)
            plt.close()

            plt.figure()
            plt.plot(band_centers, band_mean, marker='o', label="mean")
            plt.fill_between(band_centers, band_low, band_high, alpha=0.2, label=f"{int(ci*100)}% CI")
            plt.xlabel("Normalized radial frequency (band centers)")
            plt.ylabel("Fractional energy per band")
            plt.title("Spectral Energy Distribution by Bands (mean ± CI)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(out_dir / "band_energies_fractional_mean_CI.png", dpi=300)
            plt.close()

            manifest_rows.append({
                "stage": "analyze",
                "input": str(src),
                "output": str(out_dir),
                "status": "ok",
                "shape_z_y_x": f"{list(vol_zyx.shape)}"
            })
            print(f"[analyze] Outputs -> {out_dir}")

        except Exception as e:
            manifest_rows.append({
                "stage": "analyze",
                "input": str(src),
                "output": str(out_dir),
                "status": f"error: {e}",
                "shape_z_y_x": ""
            })
            print(f"[analyze][ERROR] {src}: {e}", file=sys.stderr)


# ----------------------------- CLI -----------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="FFT & Spectral Analysis for Brain MRI (NIfTI). "
                    "Mirrors input folder structure and writes clean, isolated outputs per file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--mode", choices=["fft", "analyze", "both"], required=True,
                   help="Pipeline stage to run.")
    p.add_argument("--input", required=True, type=Path,
                   help="Input root directory (spatial volumes for 'fft', FFT volumes for 'analyze').")
    p.add_argument("--output", required=True, type=Path,
                   help="Output root directory; structure is mirrored.")
    p.add_argument("--pattern", default="*.nii*",
                   help="Glob pattern for input NIfTI files.")
    p.add_argument("--n-bins", type=int, default=N_BINS_PROFILE,
                   help="Number of radial bins for the radial profile.")
    p.add_argument("--n-bands", type=int, default=N_BANDS,
                   help="Number of radial bands for band-limited energy.")
    p.add_argument("--n-boot", type=int, default=N_BOOT,
                   help="Number of bootstrap resamples.")
    p.add_argument("--ci", type=float, default=CI_LEVEL,
                   help="Confidence level for bootstrap CI (e.g., 0.95).")
    p.add_argument("--seed", type=int, default=RANDOM_SEED,
                   help="Random seed for reproducibility.")
    return p


def main():
    args = build_argparser().parse_args()

    input_root: Path = args.input.resolve()
    output_root: Path = args.output.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    manifest_rows: List[Dict] = []

    if args.mode in ("fft", "both"):
        print(f"[stage] FFT  | input={input_root} -> output={output_root}")
        stage_fft(input_root=input_root,
                  output_root=output_root,
                  pattern=args.pattern,
                  manifest_rows=manifest_rows)

        # If 'both', we should analyze the freshly written FFTs from output_root.
        if args.mode == "both":
            analyze_input = output_root
        else:
            analyze_input = None
    else:
        analyze_input = input_root

    if args.mode in ("analyze", "both"):
        print(f"[stage] ANALYZE | input={analyze_input} -> output={output_root}")
        stage_analyze(input_root=analyze_input,
                      output_root=output_root,
                      pattern=args.pattern,
                      n_bins=args.n_bins,
                      n_bands=args.n_bands,
                      n_boot=args.n_boot,
                      ci=args.ci,
                      seed=args.seed,
                      manifest_rows=manifest_rows)

    # Write manifest at the root of output
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(output_root / "manifest.csv", index=False)
    print(f"[done] Manifest saved at: {output_root / 'manifest.csv'}")


if __name__ == "__main__":
    main()
