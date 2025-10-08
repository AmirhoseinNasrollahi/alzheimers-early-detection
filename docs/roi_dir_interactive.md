# ROI Extraction CLI (Interactive, Dir-to-Dir)

Apply selected segmentation labels as a mask to registered medical images — interactively, case by case — and save the ROI-extracted outputs as `.nii.gz`.

> **What it does:**  
> - Reads registered images from `--reg` directory and segmentations from `--seg` directory  
> - **Pairs files by filename stem** (same case name, ignoring extensions)  
> - For **each matched pair**, lists available labels (with human-readable names when available, e.g., Slicer `.seg.nrrd` metadata), then **prompts you** to type label values *one by one*; type `end` to start extraction  
> - Builds a binary mask from your selected labels, **resamples the mask** to the registered image space (NearestNeighbor) if needed, **applies** it, and **saves** `<stem>-roi-extract.nii.gz` into `./roi_extraction` (configurable)  
> - Logs everything (console + `process.log`), including per-case runtime and a final summary

---

## ✨ Features

- 📁 **Directory-to-directory** processing; pairs by filename stem (`case01.nii.gz` ↔ `case01.seg.nrrd`)
- 🏷️ **Label-aware** and **interactive** per case (print label values + names; user selects labels; `end` to proceed)
- 🧠 **.seg.nrrd metadata parsing** (`SegmentN_Name`, `SegmentN_LabelValue`) to show human-readable label names
- 🔁 **Resampling** mask to registered image grid (**NearestNeighbor**) when size/spacing/origin/direction differ
- 🧪 **Robust masking** (`sitk.Mask`) preserving image meta; output is **`.nii.gz`**
- 🗂️ **Clean outputs** under `roi_extraction/` with name pattern: `<stem>-roi-extract.nii.gz`
- 🪵 **Detailed logging** (console + `process.log`), per-case timings, and a summary report
- 🧰 **Simple dependencies**: `SimpleITK`, `numpy`

---

## 🏁 Requirements

- Python 3.8+
- [SimpleITK](https://simpleitk.org/)  
- numpy

Install:

```bash
pip install SimpleITK numpy
```

(Optional) pin versions via `requirements.txt`:

```
SimpleITK>=2.3
numpy>=1.24
```

---

## 🚀 Usage

### ▶️ Run from CMD or PowerShell

#### **Example 1 – basic use (interactive):**
```bash
python roi_dir_interactive.py --reg "C:\data\regs" --seg "C:\data\segs"
```

#### **Example 2 – specify output folder and overwrite existing results:**
```bash
python roi_dir_interactive.py --reg "C:\data\regs" --seg "C:\data\segs" --outdir "C:\data\roi_extraction" --overwrite -v
```

#### **Example 3 – on Linux/macOS:**
```bash
python roi_dir_interactive.py --reg /data/regs --seg /data/segs --outdir ./roi_extraction --overwrite -v
```

---

### What you’ll see (example)

1. The tool **matches cases** by filename stem and prints how many pairs it found.  
2. For each pair, it prints a **label table** like:

```
Available labels:
  Value   Name
  0       Background
  1       Tumor_Core
  2       Edema
  4       Enhancing
```

3. You then **enter labels one by one**:
```
Label (or 'end'): 1
  ✓ added: 1 (Tumor_Core)
Label (or 'end'): 4
  ✓ added: 4 (Enhancing)
Label (or 'end'): end
```

4. It applies the mask and saves the output:
```
Saved: roi_extraction/case01-roi-extract.nii.gz (elapsed: 2.43s)
```

5. After all pairs, it prints a **summary** (processed count, successes, failures, total & average time).

---

## 🔧 Input/Output Rules

### Pairing by filename stem
- `case01.nii.gz` (registered) pairs with `case01.seg.nrrd` (segmentation), ignoring extension.
- Supported extensions: `.nii`, `.nii.gz`, `.nrrd`, `.seg.nrrd`.

### Label detection & names
- For **Slicer `.seg.nrrd`** files, label names are parsed from metadata:  
  `Segment0_Name`, `Segment0_LabelValue`, `Segment1_Name`, ...
- If names aren’t present, the script shows generic names (`Label <value>`; `Background` for `0`).

### Multi-component segmentations
- If the segmentation is **multi-component** (vector) rather than a scalar labelmap, the tool exposes:
  - `0 → Background`, `1 → Mask` (meaning “any component > 0”).
- Selecting `1` extracts “foreground” voxels; selecting other values does nothing (not a true labelmap).

### Resampling
- If mask grid differs from registered image (size/spacing/origin/direction), the mask is resampled with **NearestNeighbor** to preserve labels.

### Output naming
- Output for `case01.nii.gz` → `roi_extraction/case01-roi-extract.nii.gz`.

### Logging
- Console + file `roi_extraction/process.log`.
- Per-case timings and a final summary.

---

## ⚠️ Troubleshooting

- **“No matched (reg, seg) pairs found”**  
  Make sure filenames match by stem: `caseX.nii.gz` ↔ `caseX.seg.nrrd`.  
  Mixed casing or extra suffixes break pairing.

- **Output identical to input**  
  Your mask might be all-ones or selection included everything (including `0`).  
  Try selecting only the foreground labels (e.g., `1,4`), not `0`.

- **Different grids** (size/spacing/origin/direction)  
  Script resamples mask automatically. If results look off, verify both images are in the same physical space (registration done beforehand).

- **Huge images / memory**  
  Consider running per case; this script already processes cases sequentially.

---

## 📁 Suggested Repo Structure

```
roi-extraction-cli/
├─ roi_dir_interactive.py
├─ README.md
├─ requirements.txt
├─ LICENSE
├─ examples/
│  ├─ regs/          # sample registered images
│  └─ segs/          # matching segmentations (same stems)
└─ roi_extraction/   # output folder created automatically
```

---

## 🧾 License

MIT (recommended for research tools).

---

## 🙌 Acknowledgments

- Built with [SimpleITK](https://simpleitk.org/) and `numpy`.

---
