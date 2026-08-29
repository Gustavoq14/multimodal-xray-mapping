"""Build the reduced dataset shipped with this repository.

The raw experiment is ~4.2 GB of EDF frames (5597 detector images of 407x487
float32 pixels) and cannot be hosted on GitHub. This script extracts, from those
raw frames, the *only* information the mapping pipeline actually consumes, and
writes it to `data/` (~12 MB).

What the pipeline reads from each raw frame
-------------------------------------------
* Transmission map : a 7x7 construction kernel centred on the beam position.
* Beam-broadening  : the central 20x20 region of a 60x60 analysis window,
                     also centred on the beam position.
* Bookkeeping      : the motor coordinates (x, z) and the exposure time.

Nothing else in the 407x487 frame is used. We therefore store the 60x60 window
(a superset of both, and the window quoted in the paper) plus the motor
coordinates. Extracting a kernel of size s from the stored 60x60 window with
`extract_kernel(window, 30, 30, s)` returns exactly the same pixels as
extracting it from the full frame, so the maps are bit-identical.

Beam centre convention (important)
----------------------------------
The two code paths of the original analysis use two different reference pixels,
and the reduced dataset preserves both:

* scan frames  -> round(Center_1), round(Center_2) from the EDF header;
* direct beam  -> int() of the intensity-weighted centroid, averaged over the
                  100 direct-beam frames and computed on the FULL detector
                  frame. The centroid cannot be recomputed from a crop (the
                  asymmetric background outside the crop shifts it by ~0.9 px
                  in z), so the precomputed value is stored in metadata.json.

Usage
-----
    python tools/make_reduced_dataset.py <raw_root> [-o data]

where <raw_root> holds the two raw acquisition folders (see FOLDER_* below).
"""

import argparse
import json
import os
import re
import numpy as np

# --------------------------------------------------------------------------
# Raw acquisition folders and the file ranges used by the analysis
# --------------------------------------------------------------------------
FOLDER_SCAN = "EveLor_20260515"        # raster scans of the PHBV + BG 4% film
FOLDER_ROIS = "GusEveROIs_20260518"    # localized SAXS/WAXS at the ROIs

BEAM_SLICE = (1, 101)          # direct-beam frames used for I0 and B0
SCAN_START = 201               # first frame of the raster scan
BLOCKS = (1024, 1024)          # two interlaced scans of 1024 points each

WIN_SCAN = 60                  # beam-broadening analysis window (px)
WIN_BEAM = 80                  # direct-beam window kept for illustration (px)

SAXS_FILES = ["GusEveROIs_20260518_0_00007.edf",     # ROI1
              "GusEveROIs_20260518_0_00008.edf",     # ROI2
              "GusEveROIs_20260518_0_00009.edf"]     # ROI3
WAXS_FILE = "Resul2T.txt"      # Foxtrot azimuthal integration, one column/ROI
WAXS_COLUMNS = [3, 4, 5]       # 0-based data columns matching ROI1, ROI2, ROI3
ROI_NAMES = ["ROI1", "ROI2", "ROI3"]

# --------------------------------------------------------------------------
# Minimal EDF reader (single-image files) - avoids a fabio dependency here
# --------------------------------------------------------------------------
EDF_DTYPES = {"FloatValue": np.float32, "DoubleValue": np.float64,
              "UnsignedShort": np.uint16, "SignedInteger": np.int32,
              "UnsignedInteger": np.uint32, "SignedLong": np.int32}


def read_edf(path):
    """Return (data, header) of a single-image EDF file."""
    with open(path, "rb") as f:
        raw = f.read(8192)
        text = raw[:raw.index(b"}") + 1].decode("latin-1")
        header = {m.group(1): m.group(2).strip()
                  for m in re.finditer(r"([\w\.]+)\s*=\s*([^;]*);", text)}
        f.seek(int(header["EDF_HeaderSize"]))
        n1, n2 = int(header["Dim_1"]), int(header["Dim_2"])
        dt = np.dtype(EDF_DTYPES[header["DataType"]]).newbyteorder(
            "<" if header.get("ByteOrder", "LowByteFirst") == "LowByteFirst" else ">")
        data = np.frombuffer(f.read(n1 * n2 * dt.itemsize),
                             dtype=dt).reshape(n2, n1)
    return data, header


def centroid(image):
    """Intensity-weighted centroid (row, col) of a 2D image."""
    image = np.asarray(image, dtype=np.float64)
    total = image.sum()
    y, x = np.indices(image.shape)
    return (y * image).sum() / total, (x * image).sum() / total


def window(image, col, row, size):
    """`size`x`size` window centred on the pixel (row, col)."""
    half = size // 2
    return np.array(image[row - half:row - half + size,
                          col - half:col - half + size], dtype=np.float32)


# --------------------------------------------------------------------------
def main(raw_root, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    scan_dir = os.path.join(raw_root, FOLDER_SCAN)
    roi_dir = os.path.join(raw_root, FOLDER_ROIS)
    files = sorted(f for f in os.listdir(scan_dir) if f.endswith(".edf"))
    print(f"{len(files)} raw frames in {FOLDER_SCAN}")

    # ---------------- direct beam: centroid on the full frames --------------
    beam_files = files[BEAM_SLICE[0]:BEAM_SLICE[1]]
    cen = np.array([centroid(read_edf(os.path.join(scan_dir, f))[0])
                    for f in beam_files])                      # (N, 2) = (row, col)
    cy_mean, cx_mean = cen[:, 0].mean(), cen[:, 1].mean()
    ref_beam = (int(cx_mean), int(cy_mean))                    # (col, row)
    print(f"direct-beam centroid over {len(beam_files)} full frames: "
          f"cx={cx_mean:.4f} cy={cy_mean:.4f} -> reference pixel {ref_beam}")

    beam_stack, beam_exposure = [], []
    for f in beam_files:
        data, h = read_edf(os.path.join(scan_dir, f))
        beam_stack.append(window(data, ref_beam[0], ref_beam[1], WIN_BEAM))
        beam_exposure.append(float(h["ExposureTime"]))
    beam_stack = np.stack(beam_stack)

    np.savez_compressed(os.path.join(out_dir, "direct_beam.npz"),
                        windows=beam_stack,
                        exposure=np.array(beam_exposure, np.float32),
                        centroid_full_frame=cen.astype(np.float64))
    print(f"  direct_beam.npz  {beam_stack.shape}")

    # ---------------- raster scan: 60x60 window per scan point --------------
    n_scan = sum(BLOCKS)
    scan_files = files[SCAN_START:SCAN_START + n_scan]
    assert len(scan_files) == n_scan, "not enough scan frames"

    windows, positions, exposure = [], [], []
    ref_scan = None
    for f in scan_files:
        data, h = read_edf(os.path.join(scan_dir, f))
        col = int(round(float(h["Center_1"])))
        row = int(round(float(h["Center_2"])))
        if ref_scan is None:
            ref_scan = (col, row)
        elif (col, row) != ref_scan:
            raise ValueError(f"beam centre moves during the scan: {f}")
        windows.append(window(data, col, row, WIN_SCAN))
        positions.append((float(h["x"]), float(h["z"])))
        exposure.append(float(h["ExposureTime"]))
    windows = np.stack(windows)

    np.savez_compressed(os.path.join(out_dir, "scan_windows.npz"),
                        windows=windows,
                        positions=np.array(positions, np.float64),
                        exposure=np.array(exposure, np.float32),
                        blocks=np.array(BLOCKS, np.int32))
    print(f"  scan_windows.npz {windows.shape}  reference pixel {ref_scan}")

    # ---------------- localized SAXS frames (kept whole) --------------------
    saxs, saxs_header = [], None
    for f in SAXS_FILES:
        data, h = read_edf(os.path.join(roi_dir, f))
        saxs.append(data.astype(np.float32))
        saxs_header = h
    np.savez_compressed(os.path.join(out_dir, "saxs_rois.npz"),
                        frames=np.stack(saxs),
                        roi_names=np.array(ROI_NAMES))
    print(f"  saxs_rois.npz    {np.stack(saxs).shape}")

    # ---------------- localized WAXS profiles (already integrated) ---------
    waxs = np.loadtxt(os.path.join(roi_dir, WAXS_FILE), skiprows=1)
    two_theta = waxs[:, 0]
    cols = waxs[:, [c + 1 for c in WAXS_COLUMNS]]       # +1: column 0 is 2theta
    np.savetxt(os.path.join(out_dir, "waxs_rois.csv"),
               np.column_stack([two_theta, cols]), delimiter=",",
               header="two_theta_deg," + ",".join(ROI_NAMES), comments="")
    print(f"  waxs_rois.csv    {cols.shape[0]} points x {cols.shape[1]} ROIs")

    # ---------------- experimental parameters ------------------------------
    h0 = read_edf(os.path.join(scan_dir, scan_files[0]))[1]
    meta = {
        "description": "Reduced dataset for the multimodal X-ray mapping "
                       "notebook (PHBV + BG 4%). Derived from the raw EDF "
                       "frames by tools/make_reduced_dataset.py.",
        "instrument": {
            "name": "Xenocs Xeuss 2.0",
            "detector": h0.get("DetectorModel", ""),
            "detector_shape": [int(h0["Dim_2"]), int(h0["Dim_1"])],
            "pixel_size_m": float(h0["PSize_1"]),
            "wavelength_m": float(h0["WaveLength"]),
            "collimation": "UHR",
        },
        "scan": {
            "sample_detector_distance_m": float(h0["SampleDistance"]),
            "exposure_time_s": float(h0["ExposureTime"]),
            "motor_step_mm": 0.039,
            "blocks": list(BLOCKS),
            "n_points": int(n_scan),
            "window_px": WIN_SCAN,
            "reference_pixel_col_row": list(ref_scan),
            "reference_pixel_source": "round(Center_1), round(Center_2) from the EDF header",
            "source_folder": FOLDER_SCAN,
            "source_frames": [SCAN_START, SCAN_START + n_scan],
        },
        "direct_beam": {
            "n_frames": len(beam_files),
            "window_px": WIN_BEAM,
            "centroid_full_frame_col_row": [float(cx_mean), float(cy_mean)],
            "reference_pixel_col_row": list(ref_beam),
            "reference_pixel_source": "int() of the mean centroid over the full frames",
            "source_folder": FOLDER_SCAN,
            "source_frames": list(BEAM_SLICE),
        },
        "saxs": {
            "sample_detector_distance_m": float(saxs_header["SampleDistance"]),
            "exposure_time_s": float(saxs_header["ExposureTime"]),
            "beam_centre_col_row": [float(saxs_header["Center_1"]),
                                    float(saxs_header["Center_2"])],
            "pixel_size_m": float(saxs_header["PSize_1"]),
            "wavelength_m": float(saxs_header["WaveLength"]),
            "detector": "Pilatus200k",
            "source_folder": FOLDER_ROIS,
            "source_files": SAXS_FILES,
        },
        "waxs": {
            "source_folder": FOLDER_ROIS,
            "source_file": WAXS_FILE,
            "source_columns_zero_based": WAXS_COLUMNS,
            "note": "azimuthally integrated in Foxtrot; 2theta in degrees",
        },
        "roi_motor_coordinates_mm": {},
    }
    for name, f in zip(ROI_NAMES, SAXS_FILES):
        h = read_edf(os.path.join(roi_dir, f))[1]
        meta["roi_motor_coordinates_mm"][name] = {"x": float(h["x"]),
                                                  "z": float(h["z"])}

    with open(os.path.join(out_dir, "metadata.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)
    print("  metadata.json")

    total = sum(os.path.getsize(os.path.join(out_dir, f))
                for f in os.listdir(out_dir))
    print(f"\nreduced dataset: {total / 1e6:.1f} MB in {out_dir}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("raw_root", help="folder containing the raw acquisition folders")
    ap.add_argument("-o", "--out", default="data", help="output folder")
    args = ap.parse_args()
    main(args.raw_root, args.out)
