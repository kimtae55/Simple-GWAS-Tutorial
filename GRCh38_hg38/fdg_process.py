#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build FDG ROI matrix using AAL2 atlas (keeps roi_value=0 as Background for simple mapping).

Inputs
------
- AAL2 TSV: must contain columns 'roi_value', 'roi_name'
- AAL2 NIfTI: integer-labeled parcellation matching TSV codes
- Covariates CSV: must contain columns 'Unnamed: 0' (1-based index into volume stack) and 'id'
- FDG 4D volume stack: .npy array of shape (N_total, X, Y, Z)

Output
------
- FDG_matrix.tsv : rows=subjects, cols=ROI names (SubjectID first).
  Optionally drops the 'Background' ROI column before saving.
"""

import os
import re
import numpy as np
import pandas as pd
import nibabel as nib

# --------------------------- User config --------------------------- #
# AAL2 atlas
AAL2_TSV = "/Users/taehyo/Local_Document/data_l0ipls/FDGPET/atlas-AAL2_dseg.tsv"
AAL2_NII = "/Users/taehyo/Local_Document/data_l0ipls/FDGPET/atlas-AAL2_dseg.nii"

# FDG data
DATAPATH = "/Users/taehyo/Local_Document/data_l0ipls/FDGPET"
COV_CSV  = "ADNI_FDGPET_covariates.csv"                 # must have cols: 'Unnamed: 0', 'id'
NPY_NAME = "ADNI_FDGPET_Clinica_pet-volume_1536.npy"    # shape: (N_total, X, Y, Z)
OUT_TSV  = "FDG_matrix.tsv"

# Options
SANITIZE_NAMES = True
DROP_BACKGROUND_COLUMN = True  # drop the 'Background' column from the saved TSV
# ------------------------------------------------------------------- #

def sanitize(name: str) -> str:
    s = str(name).strip()
    s = re.sub(r"[^\w\-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

def load_aal2_labels(tsv_path: str,
                     exclude: set[int] | None = None,
                     sanitize_names: bool = True) -> tuple[list[int], list[str]]:
    """
    Load AAL2 atlas labels from TSV with columns: roi_value, roi_name.
    Keeps roi_value=0 (Background). Optionally excludes any ids in `exclude`.
    Returns two parallel lists: roi_labels (int codes) and roi_names (strings).
    """
    df = pd.read_csv(tsv_path, sep="\t")
    if not {"roi_value", "roi_name"}.issubset(df.columns):
        raise KeyError(f"{tsv_path} must have columns 'roi_value' and 'roi_name'. Got: {list(df.columns)}")

    if exclude:
        df = df[~df["roi_value"].isin(exclude)]

    # Sort by roi_value for reproducible, atlas-native ordering (0, 2001, 2002, ...)
    df = df.sort_values("roi_value")
    roi_labels = df["roi_value"].astype(int).tolist()
    roi_names = df["roi_name"].astype(str).tolist()
    if sanitize_names:
        roi_names = [sanitize(n) for n in roi_names]

    print(f"[AAL2] Loaded {len(roi_labels)} labels (including Background if present).")
    return roi_labels, roi_names

def build_index_map_from_labels(nii_path: str, roi_labels: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """
    Build a compact index map using positions in roi_labels (0..p-1), not raw label codes.
    - idx_map: same shape as NIfTI; each voxel stores position index (0..p-1)
               where position 0 corresponds to roi_labels[0] (likely Background=0).
    - counts : voxel count per position (len = p)

    This avoids huge bincount arrays due to sparse label codes like 2001, 2002, ...
    """
    vol = nib.load(nii_path)
    lbl = np.asanyarray(vol.dataobj)  # integer label volume
    p = len(roi_labels)

    # Map raw label code -> compact position (0..p-1)
    code_to_pos = {code: pos for pos, code in enumerate(roi_labels)}

    # Build index map
    idx_map = np.zeros(lbl.shape, dtype=np.int32)
    # Vectorized remap via a lookup table:
    # Because labels are sparse, do a loop over codes (few hundred) rather than fancy vectorization.
    for code, pos in code_to_pos.items():
        idx_map[lbl == code] = pos

    # Counts per position
    counts = np.bincount(idx_map.ravel(), minlength=p)
    if np.any(counts == 0):
        zero_pos = np.where(counts == 0)[0].tolist()
        zero_codes = [roi_labels[pos] for pos in zero_pos]
        raise ValueError(f"ROI(s) with zero voxels in NIfTI: positions={zero_pos}, codes={zero_codes}")

    return idx_map, counts

def compute_roi_means_subset(
    datapath: str,
    aal2_tsv: str,
    aal2_nii: str,
    cov_csv: str,
    npy_name: str,
    drop_background: bool = True,
    sanitize_names: bool = True,
    out_tsv: str | None = None
) -> pd.DataFrame:
    # 1) Load AAL2 label list (keep Background=0)
    roi_labels, roi_names = load_aal2_labels(aal2_tsv, exclude=None, sanitize_names=sanitize_names)
    p = len(roi_labels)

    # 2) Build voxel->position index map and per-ROI voxel counts
    idx_map, counts = build_index_map_from_labels(aal2_nii, roi_labels)
    flat_idx = idx_map.ravel()   # positions 0..p-1
    p_bins = p                   # bincount length

    # 3) Load subject selection (0-based indices into NPY + subject IDs)
    subjects = pd.read_csv(os.path.join(datapath, cov_csv))
    if not {"Unnamed: 0", "id"}.issubset(subjects.columns):
        raise KeyError(f"{cov_csv} must contain columns 'Unnamed: 0' and 'id'")

    subj_idx = subjects["Unnamed: 0"].to_numpy(dtype=np.int64) - 1  # convert 1-based to 0-based
    subj_ids = subjects["id"].to_numpy()

    # 4) Memory-map FDG volume stack
    vol = np.load(os.path.join(datapath, npy_name), mmap_mode="r")
    n_total = vol.shape[0]
    if (subj_idx < 0).any() or (subj_idx >= n_total).any():
        bad = subj_idx[(subj_idx < 0) | (subj_idx >= n_total)]
        raise ValueError(f"Out-of-bounds indices in 'Unnamed: 0': {bad.tolist()} (valid 1..{n_total})")
    n = subj_idx.shape[0]
    print(f"[FDG] Subjects selected: {n}")
    print(f"[FDG] Atlas ROIs (including Background): {p}")

    # 5) Compute per-ROI means (n × p) using bincount on compact positions
    X = np.empty((n, p), dtype=np.float32)
    for i, sidx in enumerate(subj_idx):
        flat_vals = vol[sidx].ravel().astype(np.float64, copy=False)
        sums = np.bincount(flat_idx, weights=flat_vals, minlength=p_bins)
        X[i, :] = (sums / counts).astype(np.float32)
        if (i + 1) % 100 == 0 or i == n - 1:
            print(f"Processed {i + 1}/{n} subjects", end="\r")
    print(f"\n[FDG] ROI matrix shape: {X.shape} (n × p)")

    # 6) Build DataFrame; optionally drop Background column
    df = pd.DataFrame(X, columns=roi_names)
    df.insert(0, "SubjectID", subj_ids)

    if drop_background:
        # Background is the row where roi_labels == 0; find its position
        if 0 in roi_labels:
            bg_pos = roi_labels.index(0)
            bg_name = roi_names[bg_pos]
            if bg_name in df.columns:
                df = df.drop(columns=[bg_name])
                print(f"[FDG] Dropped Background column: '{bg_name}'")
        else:
            print("[FDG] Background (0) not found in roi_labels; nothing dropped.")

    # 7) Save
    if out_tsv is not None:
        out_path = os.path.join(datapath, out_tsv)
        df.to_csv(out_path, sep="\t", index=False, na_rep="NA")
        print(f"[FDG] Saved ROI matrix TSV: {out_path}")

    return df

if __name__ == "__main__":
    df = compute_roi_means_subset(
        datapath=DATAPATH,
        aal2_tsv=AAL2_TSV,
        aal2_nii=AAL2_NII,
        cov_csv=COV_CSV,
        npy_name=NPY_NAME,
        drop_background=DROP_BACKGROUND_COLUMN,
        sanitize_names=SANITIZE_NAMES,
        out_tsv=OUT_TSV
    )
    print(df.iloc[:5, :6])  # preview