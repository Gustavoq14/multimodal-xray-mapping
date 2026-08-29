# Multimodal X-ray mapping from the transmitted beam

Reference implementation of the mapping workflow described in

> **Multimodal X-ray mapping from the transmitted beam to guide localized SAXS/WAXS on a
> laboratory instrument.**
> G. Quintero-Marquina, F. Britto, E. Pacheco, B. Araoz, D. G. Lamas, C. Huck-Iriart.

The method turns a commercial laboratory SAXS/WAXS instrument (Xenocs Xeuss 2.0) into a
two-dimensional mapping system without any instrumental modification. It uses the intrinsic
divergence of the beam as a geometric magnifier, and recombines two interlaced raster scans
with the *drizzle* algorithm to reach a sampling step of 19.5 µm. From a single acquisition
the same transmitted-beam images yield three complementary maps — transmission, absorption
and beam broadening — the last of which is used to locate heterogeneities and to place the
regions of interest (ROIs) at which localized SAXS/WAXS is subsequently measured.

The notebook covers the **PHBV + Bioglass 4 wt %** system.

## Contents

```
multimodal_mapping_PHBV_BG4.ipynb   the workflow, from the scan data to the SAXS/WAXS profiles
data/                               reduced dataset (~11 MB), enough to run the notebook
tools/make_reduced_dataset.py       how data/ was derived from the raw detector frames
requirements.txt                    Python dependencies
```

## Running it

```bash
pip install -r requirements.txt
jupyter lab multimodal_mapping_PHBV_BG4.ipynb
```

Everything the notebook needs is in this repository; no access to the raw experimental data
is required. It runs end to end in a few seconds.

## What the notebook does

```
raster scans  ->  construction kernels  ->  transmission and beam-broadening maps
              ->  drizzle (sub-pixel recombination)
              ->  transmission T, absorption mu*t, beam broadening B
              ->  ROI selection on B  ->  localized SAXS / WAXS profiles
```

| map | symbol | definition |
|---|---|---|
| transmission | *T* | *I* / *I*₀ |
| absorption | *μt* | −ln *T* |
| beam broadening | *B* | (σ/*I*_c) / (σ₀/*I*_c,₀) |

## The data

The raw experiment is about **4.2 GB**: 5593 detector frames of 407 × 487 float32 pixels
for the raster scans, plus the localized SAXS/WAXS acquisitions. That is impractical to
host here, so this repository ships a **reduced dataset of ~11 MB** in `data/`.

The reduction is not arbitrary. Of each 407 × 487 detector frame the pipeline only ever
reads a small window centred on the beam: a 7 × 7 construction kernel for the transmission
map, and the central 20 × 20 region of the analysis window for the beam-broadening map.
`data/` therefore stores the 60 × 60 window quoted in the article — a superset of both —
together with the motor coordinates and exposure times. Cropping a kernel from the stored
window is the same arithmetic as cropping it from the full frame, so **the maps built from
the reduced dataset are bit-identical to those built from the raw frames** (verified
against the full dataset).

| file | contents |
|---|---|
| `scan_windows.npz` | `windows` (2048, 60, 60) detector windows, one per scan point; `positions` (2048, 2) motor coordinates (x, z) in mm; `exposure` (2048,) in s; `blocks` the two interlaced scans |
| `direct_beam.npz` | `windows` (100, 80, 80) direct-beam windows; `exposure`; `centroid_full_frame` (100, 2) per-frame centroids computed on the full frames |
| `saxs_rois.npz` | `frames` (3, 407, 487) 2D SAXS patterns measured at ROI1–ROI3 |
| `waxs_rois.csv` | azimuthally integrated WAXS profiles at ROI1–ROI3 (2θ in degrees) |
| `metadata.json` | instrument, scan, direct-beam, SAXS and WAXS parameters, and the measured ROI motor coordinates |

### One step that cannot be recomputed from the reduced data

The beam centroid used to place the direct-beam kernel is an intensity-weighted centroid of
the **full** detector frame. Recomputing it on an 80 × 80 crop shifts it by ≈ 0.9 px in *z*,
because the asymmetric background outside the crop contributes to the full-frame value. The
precomputed per-frame centroids are therefore stored in `direct_beam.npz`, and the windows
were cropped around the corresponding reference pixel so that every kernel extraction
downstream is exact. The notebook shows this comparison explicitly.

### Regenerating `data/` from the raw frames

`tools/make_reduced_dataset.py` documents the exact provenance of every array (source
folder, frame range, reference pixel and the convention used for it). Pointing it at the
raw acquisition folders regenerates `data/` from scratch:

```bash
python tools/make_reduced_dataset.py /path/to/raw_data -o data
```

The raw frames are kept outside this repository and are available from the authors.

## The method, step by step

The animations below illustrate the three computational steps of the method. They were
produced with the animation routines of the full analysis notebook (not included here,
since they are not needed to reproduce the maps).

> **Note.** These animations show the **PLA + SBA-15 4 wt %** sample, the other system
> studied in the article, whose larger scan (2025 points per raster) makes the mechanism
> easier to follow. The notebook in this repository processes the PHBV + BG 4 % sample,
> so the maps it produces look different — the procedure is identical.

**1. Building the transmission map, kernel by kernel.** Each scan point contributes a
7 × 7 construction kernel, deposited on the common grid at the motor position where it was
acquired. Neighbouring kernels overlap, so each map pixel accumulates several
contributions and is divided by its coverage *N*(**p**) (right panel). The map is finally
cropped to the region where *N* = 49, the only one with complete coverage.

![Transmission map construction](assets/Transmission_map.gif)

**2. Building the beam-broadening map, point by point.** Here each scan point yields a
single scalar instead of a kernel: the radial second moment σ of the transmitted beam
divided by the intensity at its centroid *I*<sub>c</sub>, both measured on the analysis
window (right panel, with σ drawn as the dashed circle).

![Beam-broadening map construction](assets/Beam_broadening_map.gif)

**3. Sub-pixel recombination by drizzle.** The two interlaced scans, offset by half a motor
step, are deposited on a grid of half the step. Each input pixel is shrunk to `pixfrac`
= 0.6 of its size before being dropped, which is what turns the half-step offset into
distinct output positions instead of an average. The right panel shows the accumulated
weight converging to a uniform *w* = 0.72.

![Drizzle recombination](assets/Drizzle_map.gif)

## Citation

If you use this code, please cite the article above.
