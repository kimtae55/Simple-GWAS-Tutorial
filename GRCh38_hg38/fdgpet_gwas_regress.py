import pandas as pd
import numpy as np
import re

# ===== Paths =====
fam_path  = "/Users/taehyo/Local_Document/data_l0ipls/ADNI_MERGED_FINAL/ADNI_qc_final.fam"
adni_path = "/Users/taehyo/Local_Document/data_l0ipls/Demographic/ADNIMERGE.csv"
fdg_path  = "/Users/taehyo/Local_Document/data_l0ipls/FDGPET/ADNI_FDGPET_covariates.csv"

# ===== Helpers =====
def str_trim_upper(s):
    if pd.isna(s): return s
    return str(s).strip().upper()

def fdg_id_to_ptid(id_str):
    if pd.isna(id_str): return np.nan
    m = re.match(r"^sub-ADNI(\d{3})S(\d+)$", str(id_str))
    return f"{m.group(1)}_S_{m.group(2)}".upper() if m else str(id_str).upper()

def summarize_cat(df, col):
    vc = df[col].value_counts(dropna=False)
    total = int(vc.sum())
    out = (pd.DataFrame({"level": vc.index, "N": vc.values})
             .assign(pct=lambda x: (100 * x["N"] / total).round(1).astype(str) + "%"))
    return out.sort_values("N", ascending=False, ignore_index=True)

# ===== 1) Load fam and normalize to PTID_like =====
fam = pd.read_csv(
    fam_path, sep=r"\s+", header=None,
    names=["FID","IID","MID","PID","SEX","PHENO"], engine="python"
)
fam["PTID_like"] = fam["IID"].astype(str).str.replace(r"^[^_]+_", "", regex=True).str.upper()

# ===== 2) Load ADNIMERGE (baseline only) =====
ad = pd.read_csv(adni_path, usecols=["PTID","VISCODE","COLPROT","DX_bl","AGE","PTGENDER","PTEDUCAT","PTETHCAT","PTRACCAT","PTMARRY"])
ad["PTID"] = ad["PTID"].apply(str_trim_upper)
ad = ad[ad["VISCODE"] == "bl"].drop(columns=["VISCODE"]).set_index("PTID")

# ===== 3) Load FDG covariates, map to PTID, keep original FDG id =====
fdg = pd.read_csv(fdg_path)
fdg_map = (fdg.assign(PTID=lambda d: d["id"].apply(fdg_id_to_ptid))
              .drop_duplicates("PTID")[["PTID","id"]]
              .rename(columns={"id":"FDG_ID"})
              .set_index("PTID"))

# ===== 4) Keep subjects present in BOTH ADNIMERGE and FDG =====
keep = fam["PTID_like"].isin(ad.index.intersection(fdg_map.index))
fam_both = fam.loc[keep, ["FID","IID","PTID_like","SEX"]].copy()

# ===== 5) Join by PTID_like (no duplicate PTID columns) =====
m = (fam_both
     .set_index("PTID_like")
     .join(ad, how="left")
     .join(fdg_map, how="left")
     .reset_index())

# ===== 6) Clean types =====
m["AGE"] = pd.to_numeric(m["AGE"], errors="coerce")
m["PTEDUCAT"] = pd.to_numeric(m["PTEDUCAT"], errors="coerce")
for c in ["COLPROT","PTGENDER","PTETHCAT","PTRACCAT","PTMARRY","DX_bl"]:
    m[c] = m[c].astype(str).str.strip().replace({"nan": np.nan})

# ===== 7) Final covariate table =====
covar_final = m.rename(columns={"PTID_like":"PTID"})[
    ["FID","IID","FDG_ID","SEX","AGE","PTEDUCAT",
     "COLPROT","PTGENDER","PTETHCAT","PTRACCAT","PTMARRY","DX_bl"]
].copy()

# ===== 8) Save + quick distributions =====
outfile = "covar_ADNI_baseline_fdgpet_gwas_adnimerge.tsv"
covar_final.to_csv(outfile, sep="\t", index=False, na_rep="NA", quoting=3)

print("\n== Final distributions (based on rows saved to file) ==")
print("\n-- DX_bl --");    print(summarize_cat(covar_final, "DX_bl"))
print("\n-- COLPROT --"); print(summarize_cat(covar_final, "COLPROT"))
print("\n-- PTGENDER --");print(summarize_cat(covar_final, "PTGENDER"))
print(f"\nSaved: {outfile}")


######################################################################################################## 
# Residualization of FDG_matrix.tsv and GWAS_matrix.tsv
# output: (n_cn, p) and (n_ad, p) data matrix for SNP, and (n_cn, q) and (n_ad, q) matrix for FDG PET 
#         All datasets are column-centered with mean 0, and scaled to have unit variance 
######################################################################################################## 
# Matrix variables
''' 
# FDG_matrix.tsv
SubjectID
sub-ADNI014S2185
sub-ADNI035S4082
sub-ADNI072S4131
sub-ADNI062S0768
sub-ADNI094S4434
sub-ADNI094S0489
sub-ADNI130S5059
sub-ADNI013S4579
sub-ADNI011S0861

# GWAS_matrix.tsv
FID	IID
1_014_S_0520	1_014_S_0520
1_137_S_4303	1_137_S_4303
2_002_S_4746	2_002_S_4746
2_002_S_4799	2_002_S_4799
2_002_S_5018	2_002_S_5018
2_002_S_5178	2_002_S_5178
2_002_S_5256	2_002_S_5256
2_003_S_4524	2_003_S_4524
2_003_S_4644	2_003_S_4644
2_003_S_4872	2_003_S_4872
2_003_S_4892	2_003_S_4892

# merged information based on ADNIMERGE
FID	IID	FDG_ID	SEX	AGE	PTEDUCAT	COLPROT	PTGENDER	PTETHCAT	PTRACCAT	PTMARRY	DX_bl
1_137_S_4303	1_137_S_4303	sub-ADNI137S4303	2	80.2	20	ADNI2	Female	Not Hisp/Latino	White	Widowed	LMCI
2_002_S_4746	2_002_S_4746	sub-ADNI002S4746	2	71.2	16	ADNI2	Female	Not Hisp/Latino	White	Married	LMCI
'''
########################################################################################################
# ===== Imports =====
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

# ===== Paths =====
FDG_TSV   = "FDG_matrix.tsv"
GWAS_TSV  = "GWAS_matrix.tsv"
COVAR_TSV = "covar_ADNI_baseline_fdgpet_gwas_adnimerge.tsv"
OUTDIR    = "SCCA_analysis"
os.makedirs(OUTDIR, exist_ok=True)

def to_ptid_like_from_sub_adni(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.upper()
    return s.str.replace(r"^SUB-ADNI(\d{3})S(\d+)$", r"\1_S_\2", regex=True)

def to_ptid_like_from_iid(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip().str.upper()
    s = s.str.replace(r"^\d+_", "", regex=True)
    s = s.str.replace(r"^(\d{3})S(\d+)$", r"\1_S_\2", regex=True)
    s = s.str.replace(r"^(\d{3})[_-]?S[_-]?(\d+)$", r"\1_S_\2", regex=True)
    return s

def residualize(data: pd.DataFrame, covariates: pd.DataFrame) -> pd.DataFrame:
    X = covariates.to_numpy()
    Y = data.to_numpy()
    XtX = X.T @ X
    XtX_pinv = np.linalg.pinv(XtX)
    H = X @ XtX_pinv @ X.T
    P = np.eye(H.shape[0]) - H
    R = P @ Y
    return pd.DataFrame(R, index=data.index, columns=data.columns)

def zscore_within_group(D: pd.DataFrame) -> pd.DataFrame:
    mu = D.mean(axis=0)
    sd = D.std(axis=0, ddof=1).replace(0, 1.0)
    return (D - mu) / sd

FDG_raw  = pd.read_csv(FDG_TSV, sep="\t")
GWAS_raw = pd.read_csv(GWAS_TSV, sep="\t")
COV      = pd.read_csv(COVAR_TSV, sep="\t")

FDG_raw["PTID_like"]  = to_ptid_like_from_sub_adni(FDG_raw["SubjectID"])
GWAS_raw["PTID_like"] = to_ptid_like_from_iid(GWAS_raw["IID"])
COV["PTID_like"]      = to_ptid_like_from_iid(COV["IID"])

common = set(FDG_raw["PTID_like"]) & set(GWAS_raw["PTID_like"]) & set(COV["PTID_like"])

FDG = (FDG_raw[FDG_raw["PTID_like"].isin(common)]
       .drop_duplicates("PTID_like").set_index("PTID_like").drop(columns=["SubjectID"]))
GWAS = (GWAS_raw[GWAS_raw["PTID_like"].isin(common)]
        .drop_duplicates("PTID_like").set_index("PTID_like").drop(columns=["FID","IID"]))
COV = (COV[COV["PTID_like"].isin(common)]
       .drop_duplicates("PTID_like").set_index("PTID_like"))

FDG  = FDG.apply(pd.to_numeric, errors="raise")
GWAS = GWAS.apply(pd.to_numeric, errors="raise")

cov_keep = ["COLPROT","AGE","PTEDUCAT","PTRACCAT","PTETHCAT","PTMARRY","PTGENDER","DX_bl"]
C = COV[cov_keep].copy()
C["AGE"] = pd.to_numeric(C["AGE"], errors="raise")
C["PTEDUCAT"] = pd.to_numeric(C["PTEDUCAT"], errors="raise")

X = pd.get_dummies(
    C[["COLPROT","PTRACCAT","PTETHCAT","PTMARRY","PTGENDER","AGE","PTEDUCAT"]],
    drop_first=True
).astype(float)
X = sm.add_constant(X, has_constant="add")

FDG  = FDG.loc[X.index].astype(float)
GWAS = GWAS.loc[X.index].astype(float)

FDG_res  = residualize(FDG,  X)
GWAS_res = residualize(GWAS, X)

dx = C["DX_bl"].astype(str).str.upper().str.strip()
is_cn = dx.eq("CN")
is_ad = dx.eq("AD")

FDG_cn,  FDG_ad  = FDG_res.loc[is_cn.values],  FDG_res.loc[is_ad.values]
GWAS_cn, GWAS_ad = GWAS_res.loc[is_cn.values], GWAS_res.loc[is_ad.values]

FDG_cn_z  = zscore_within_group(FDG_cn)
FDG_ad_z  = zscore_within_group(FDG_ad)
GWAS_cn_z = zscore_within_group(GWAS_cn)
GWAS_ad_z = zscore_within_group(GWAS_ad)

FDG_cn_z.to_csv(os.path.join(OUTDIR, "fdg_cn.tsv"),  sep="\t")
FDG_ad_z.to_csv(os.path.join(OUTDIR, "fdg_ad.tsv"),  sep="\t")
GWAS_cn_z.to_csv(os.path.join(OUTDIR, "gwas_cn.tsv"), sep="\t")
GWAS_ad_z.to_csv(os.path.join(OUTDIR, "gwas_ad.tsv"), sep="\t")

print("\n== Dataset shapes ==")
print("FDG_cn_z:", FDG_cn_z.shape)
print("FDG_ad_z:", FDG_ad_z.shape)
print("GWAS_cn_z:", GWAS_cn_z.shape)
print("GWAS_ad_z:", GWAS_ad_z.shape)