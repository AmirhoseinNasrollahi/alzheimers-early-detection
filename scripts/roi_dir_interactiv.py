#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Interactive ROI Extraction
===================================================

Reads all registered images from a directory and matches each one with a
segmentation (same filename stem) from another directory. For each matched case,
it *interactively* lists available labels (with human-readable names if present,
e.g., Slicer .seg.nrrd metadata), asks the user to enter label values one by one
(type 'end' to finish), builds a mask from the selected labels, applies it to the
registered image, and saves the masked image as <stem>-roi-extract.nii.gz under
the output directory (default: ./roi_extraction).

Key features
------------
- Pairing by filename stem (case-insensitive), ignoring dual extensions (.nii.gz, .seg.nrrd, etc.)
- .seg.nrrd label-name parsing (SegmentN_Name / SegmentN_LabelValue)
- Binary mask creation from selected integer labels
- Automatic resampling of the mask to the registered image grid (NearestNeighbor)
- Detailed logging (console + file), per-case timing, and a final summary

Dependencies
------------
- numpy
- SimpleITK

Usage
-----
    python roi_dir_interactive.py --reg /path/to/regs --seg /path/to/segs
    # optional:
    python roi_dir_interactive.py --reg regs --seg segs --outdir results/roi_extraction --overwrite -v

Notes
-----
- Interactive per-case input: enter label integers one by one; type 'end' to begin processing.
- If a segmentation is multi-component (e.g., vector image), this script cannot reliably map per-label values;
  it will treat any >0 as foreground only if you select label '1' (shown as 'Mask'). For typical labelmaps
  (scalar), it lists actual integer label values present.

"""

import argparse
import sys
import time
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import SimpleITK as sitk
import re


SUPPORTED_IMG_EXT = {".nii", ".nii.gz", ".nrrd", ".seg.nrrd"}


# ------------------------------ Logging ------------------------------

def setup_logger(outdir: Path, verbose: bool) -> logging.Logger:
    """Configure console + file logger."""
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "process.log"

    logger = logging.getLogger("roi_extraction")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(ch)

    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    logger.debug(f"Logger initialized. Log file: {log_file}")
    return logger


# ------------------------------ File utils ------------------------------

def is_image_file(p: Path) -> bool:
    name = p.name.lower()
    return any(name.endswith(ext) for ext in SUPPORTED_IMG_EXT)

def list_images(root: Path) -> List[Path]:
    """List single level only (no recursion) — directory must exist."""
    if not root.exists():
        return []
    if root.is_file():
        return [root] if is_image_file(root) else []
    return [p for p in root.iterdir() if p.is_file() and is_image_file(p)]

def stem_key(path: Path) -> str:
    """
    Pairing key from filename, ignoring dual extensions.
    E.g., 'case01.nii.gz' -> 'case01', 'brain.seg.nrrd' -> 'brain'.
    """
    name = path.name.lower()
    for suffix in (".nii.gz", ".seg.nrrd", ".nrrd", ".nii"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem.lower()

def pair_by_stem(reg_paths: List[Path], seg_paths: List[Path]) -> List[Tuple[Path, Path]]:
    """Return list of (reg, seg) pairs where stems match (case-insensitive)."""
    reg_map = {stem_key(rp): rp for rp in reg_paths}
    seg_map = {stem_key(sp): sp for sp in seg_paths}
    keys = sorted(set(reg_map.keys()) & set(seg_map.keys()))
    return [(reg_map[k], seg_map[k]) for k in keys]


# ------------------------------ I/O helpers ------------------------------

def sitk_read(path: Path) -> sitk.Image:
    return sitk.ReadImage(str(path))

def sitk_write(img: sitk.Image, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(img, str(path))


# ------------------------------ .seg.nrrd label metadata ------------------------------

def parse_seg_labels_from_metadata(seg_img: sitk.Image, logger: logging.Logger) -> Dict[int, str]:
    """
    Extract label names from Slicer .seg.nrrd metadata:
        Segment0_Name, Segment0_LabelValue, Segment1_Name, Segment1_LabelValue, ...
    Returns dict {LabelValue(int): SegmentName(str)}.
    """
    mapping: Dict[int, str] = {}
    keys = seg_img.GetMetaDataKeys()
    segments: Dict[str, Dict[str, str]] = {}
    for k in keys:
        if "Segment" in k:
            m = re.match(r"(Segment\d+)_(Name|LabelValue)", k)
            if m:
                seg_id = m.group(1)
                field = m.group(2)
                val = seg_img.GetMetaData(k)
                segments.setdefault(seg_id, {})[field] = val

    for seg_id, fields in segments.items():
        if "LabelValue" in fields:
            try:
                lv = int(fields["LabelValue"])
            except ValueError:
                continue
            name = fields.get("Name", f"Label {lv}")
            mapping[lv] = name

    if mapping:
        logger.debug(f"Parsed {len(mapping)} label(s) from .seg.nrrd metadata.")
    else:
        logger.debug("No .seg.nrrd label metadata found.")
    return mapping


# ------------------------------ Label discovery & mask making ------------------------------

def discover_labels(seg_img: sitk.Image, logger: logging.Logger) -> Tuple[List[int], Dict[int, str], bool]:
    """
    Returns (sorted_label_values, label_name_map, is_multicomponent).
    - If multicomponent (e.g., vector mask), it isn't a classic integer labelmap.
      In that case we expose {0,1} with names {Background, Mask}.
    """
    comps = seg_img.GetNumberOfComponentsPerPixel()
    arr = sitk.GetArrayFromImage(seg_img)

    if comps and comps > 1:
        # Multi-component: treat any component > 0 as 'Mask'
        labels = sorted(np.unique((arr > 0).any(axis=-1).astype(np.uint8)).tolist())
        name_map = {0: "Background", 1: "Mask"}
        return labels, name_map, True

    # Scalar (classic labelmap)
    label_name_map = parse_seg_labels_from_metadata(seg_img, logger)
    unique_vals = np.unique(arr).astype(int)
    labels = sorted(unique_vals.tolist())

    # Ensure we have a name for each value
    for lv in labels:
        if lv not in label_name_map:
            label_name_map[lv] = "Background" if lv == 0 else f"Label {lv}"

    return labels, label_name_map, False

def build_mask_from_selected_labels(seg_img: sitk.Image, selected: Set[int], is_multicomponent: bool) -> sitk.Image:
    """
    Build a binary mask (0/1) from selected label values.
    - For multi-component seg: selecting {1} yields any-component>0 foreground; otherwise zero mask.
    - For scalar labelmap: mask = 1 where seg in selected labels.
    """
    arr = sitk.GetArrayFromImage(seg_img)
    if is_multicomponent:
        bin_arr = (arr > 0).any(axis=-1).astype(np.uint8)
        if 1 not in selected:
            bin_arr[:] = 0
    else:
        sel = np.isin(arr, list(selected))
        bin_arr = sel.astype(np.uint8)

    out = sitk.GetImageFromArray(bin_arr)
    out.CopyInformation(seg_img)
    return out


# ------------------------------ Resampling & masking ------------------------------

def resample_to_reference(moving: sitk.Image, reference: sitk.Image, is_label: bool = True) -> sitk.Image:
    """Resample moving image to reference grid (NearestNeighbor for labels)."""
    interpolator = sitk.sitkNearestNeighbor if is_label else sitk.sitkLinear
    return sitk.Resample(
        moving,
        reference,
        sitk.Transform(),
        interpolator,
        reference.GetOrigin(),
        reference.GetSpacing(),
        reference.GetDirection(),
        0,
        moving.GetPixelID(),
    )

def apply_mask_to_image(image: sitk.Image, mask: sitk.Image) -> sitk.Image:
    """Keep intensities where mask != 0, else set to 0 (preserve meta)."""
    out = sitk.Mask(image, sitk.Cast(mask, sitk.sitkUInt8))
    out.CopyInformation(image)
    return out


# ------------------------------ Per-case interactive workflow ------------------------------

def prompt_labels_interactively(labels: List[int], name_map: Dict[int, str]) -> Set[int]:
    """
    Show available labels and prompt the user to type integers one by one.
    Type 'end' to finish. Returns the selected label set (can be empty).
    """
    print("\nAvailable labels:")
    print("  Value\tName")
    for lv in labels:
        print(f"  {lv}\t{name_map.get(lv, f'Label {lv}')}")
    print("\nEnter labels to extract (one per line). Type 'end' to start extraction.")
    selected: Set[int] = set()
    while True:
        raw = input("Label (or 'end'): ").strip()
        if raw.lower() == "end":
            break
        if raw.startswith("-") and raw[1:].isdigit() or raw.isdigit():
            v = int(raw)
            if v in labels:
                selected.add(v)
                print(f"  ✓ added: {v} ({name_map.get(v, f'Label {v}')})")
            else:
                print(f"  ! label {v} not found in this segmentation.")
        else:
            print("  ! please enter an integer label or 'end'.")
    return selected


# ------------------------------ Case processing ------------------------------

def process_case(reg_path: Path,
                 seg_path: Path,
                 outdir: Path,
                 overwrite: bool,
                 logger: logging.Logger) -> Tuple[Optional[Path], float, Optional[str]]:
    """
    Process a single (reg, seg) pair in interactive mode.
    Returns (output_path_or_None, elapsed_seconds, error_message_or_None).
    """
    t0 = time.time()
    try:
        logger.info(f"Processing case:\n  REG: {reg_path}\n  SEG: {seg_path}")
        reg = sitk_read(reg_path)
        seg = sitk_read(seg_path)

        # Discover labels and (optionally) names
        labels, name_map, is_multi = discover_labels(seg, logger)

        # Ask user which labels to extract
        selected = prompt_labels_interactively(labels, name_map)
        if not selected:
            msg = "No labels selected; skipping."
            logger.warning(msg)
            return None, time.time() - t0, msg

        # Build binary mask
        mask = build_mask_from_selected_labels(seg, selected, is_multi)

        # Align mask to reg space if needed
        needs_resample = (
            reg.GetSize() != mask.GetSize()
            or not np.allclose(reg.GetSpacing(), mask.GetSpacing())
            or reg.GetDirection() != mask.GetDirection()
            or not np.allclose(reg.GetOrigin(), mask.GetOrigin())
        )
        if needs_resample:
            logger.debug("Resampling mask to registered image space (NearestNeighbor).")
            mask = resample_to_reference(mask, reg, is_label=True)
        else:
            logger.debug("Mask already aligned with registered image space.")

        # Apply mask
        masked = apply_mask_to_image(reg, mask)

        # Save output as <stem>-roi-extract.nii.gz
        outdir.mkdir(parents=True, exist_ok=True)
        out_name = f"{stem_key(reg_path)}-roi-extract.nii.gz"
        out_path = outdir / out_name

        if out_path.exists() and not overwrite:
            msg = f"Output exists; use --overwrite to replace: {out_path}"
            logger.warning(msg)
            return None, time.time() - t0, msg

        sitk_write(masked, out_path)

        # Sanity note
        if np.array_equal(sitk.GetArrayFromImage(reg), sitk.GetArrayFromImage(masked)):
            logger.warning("Masked image equals original (mask may be all-ones).")

        elapsed = time.time() - t0
        logger.info(f"Saved: {out_path} (elapsed: {elapsed:.2f}s)")
        return out_path, elapsed, None

    except Exception as e:
        elapsed = time.time() - t0
        logger.exception(f"Error while processing: {e}")
        return None, elapsed, str(e)


# ------------------------------ CLI ------------------------------

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Interactive directory-to-directory ROI extractor. "
                    "Pairs registered images with segmentations by filename stem and prompts "
                    "for labels to extract for each case.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--reg", required=True, type=Path,
                   help="Directory containing registered images (or a single file).")
    p.add_argument("--seg", required=True, type=Path,
                   help="Directory containing segmentations (or a single file).")
    p.add_argument("--outdir", type=Path, default=Path("./roi_extraction"),
                   help="Output directory for <stem>-roi-extract.nii.gz and logs.")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing outputs if present.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose console logging.")
    return p

def main():
    args = build_argparser().parse_args()
    logger = setup_logger(args.outdir, args.verbose)

    # Gather inputs (single-file also allowed; but tool is designed for dirs)
    reg_items = list_images(args.reg)
    seg_items = list_images(args.seg)

    if not reg_items:
        logger.error(f"No images found in --reg: {args.reg}")
        sys.exit(1)
    if not seg_items:
        logger.error(f"No images found in --seg: {args.seg}")
        sys.exit(1)

    pairs = pair_by_stem(reg_items, seg_items)
    if not pairs:
        logger.error("No matched (reg, seg) pairs found (by filename stem).")
        logger.info(f"Found REG candidates: {len(reg_items)}, SEG candidates: {len(seg_items)}")
        sys.exit(1)

    logger.info(f"Matched pairs: {len(pairs)}")
    for idx, (rp, sp) in enumerate(pairs, 1):
        logger.debug(f"PAIR {idx}: {rp.name}  <->  {sp.name}")

    # Process all pairs interactively
    results = []
    t_all0 = time.time()
    for rp, sp in pairs:
        out_path, elapsed, err = process_case(rp, sp, args.outdir, args.overwrite, logger)
        results.append((out_path, elapsed, err))

    # Summary
    total_elapsed = time.time() - t_all0
    n_total = len(results)
    n_ok = sum(1 for r in results if r[0] is not None and r[2] is None)
    n_err = sum(1 for r in results if r[2] is not None)
    avg_time = (sum(r[1] for r in results) / n_total) if n_total else 0.0

    logger.info("---- SUMMARY ----")
    logger.info(f"Processed: {n_total} case(s)")
    logger.info(f"Succeeded: {n_ok}")
    logger.info(f"Failed:    {n_err}")
    logger.info(f"Total time: {total_elapsed:.2f}s | Avg per case: {avg_time:.2f}s")

    if n_err:
        logger.info("Failed cases:")
        for (out_path, elapsed, err), (rp, sp) in zip(results, pairs):
            if err:
                logger.info(f"- REG={rp.name}, SEG={sp.name}, elapsed={elapsed:.2f}s, error={err}")


if __name__ == "__main__":
    main()
