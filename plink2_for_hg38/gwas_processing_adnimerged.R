suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ========= Paths (edit if needed) =========
fam_path  <- "/Users/taehyo/Local_Document/data_l0ipls/ADNI_MERGED_FINAL/ADNI_qc_final.fam"
adnm_path <- "/Users/taehyo/Local_Document/data_l0ipls/Demographic/ADNIMERGE.csv"
evec_path <- "/Users/taehyo/Local_Document/data_l0ipls/ADNI_MERGED_FINAL/ADNI_pca.eigenvec"  # plink2 --pca 10

# ========= Small helpers =========
to_chr     <- function(x){ x <- as.character(x); x[is.na(x)] <- NA_character_; trimws(x) }
norm_lower <- function(x){ tolower(to_chr(x)) }

# ========= 1) Load FAM and derive PTID-like key =========
fam <- fread(fam_path, header = FALSE,
             col.names = c("FID","IID","MID","PID","SEX_FAM","PHENO"))
fam[, PTID_like := toupper(sub("^[^_]+_", "", IID))]  # "1_014_S_0520" → "014_S_0520"

# ========= 2) Load ADNIMERGE baseline rows (keep needed columns) =========
need_cols <- c("PTID","VISCODE","COLPROT","EXAMDATE","DX_bl","AGE",
               "PTGENDER","PTEDUCAT","PTETHCAT","PTRACCAT","PTMARRY")
ad <- fread(adnm_path, select = need_cols)

# Normalize keys/strings
ad[, PTID := toupper(str_trim(PTID))]
ad[, VISCODE := norm_lower(VISCODE)]
ad <- ad[VISCODE == "bl"]  # baseline only

# If multiple baseline-like rows exist, keep the earliest EXAMDATE
setorder(ad, PTID, EXAMDATE)
ad <- ad[!duplicated(PTID)]

# ========= 3) Anchor: IDs present in BOTH fam and ADNIMERGE(bl) =========
anchor <- merge(
  fam[, .(FID, IID, PTID_like)],
  ad,
  by.x = "PTID_like", by.y = "PTID",
  all = FALSE
)

# ========= 4) Phenotype (bind as a COLUMN to avoid order bugs) =========
anchor[, DX_bl_s := norm_lower(DX_bl)]
anchor[, DIAG01 := fcase(
  DX_bl_s == "cn", 1L,
  DX_bl_s == "ad", 2L,
  default = -9L
)]
anchor[, DX_bl_s := NULL]

cat("\n== Check DX_bl → DIAG01 mapping (pre-merge) ==\n")
print(anchor[, .N, by = .(DX_bl, DIAG01)][order(DX_bl, DIAG01)])

# Carry DIAG01 forward
m <- copy(anchor)

# ========= 5) Covariates (from ADNIMERGE) =========
# AGE & PTEDUCAT numeric
m[, AGE      := suppressWarnings(as.numeric(AGE))]
m[, PTEDUCAT := suppressWarnings(as.numeric(PTEDUCAT))]

# Sex (from ADNIMERGE PTGENDER strings): Female=1, Male=0, else NA
pg <- norm_lower(m$PTGENDER)
m[, PTGENDER_Female := fifelse(pg == "female", 1L,
                               fifelse(pg == "male", 0L, NA_integer_))]

# Ethnicity (PTETHCAT strings): Hisp/Latino=1, Not Hisp/Latino=0, else NA
pe <- norm_lower(m$PTETHCAT)
m[, PTETHCAT_HispLatino := fifelse(pe == "hisp/latino", 1L,
                                   fifelse(pe == "not hisp/latino", 0L, NA_integer_))]

# Race (PTRACCAT strings; White is reference → no dummy)
pr <- norm_lower(m$PTRACCAT)
m[, `:=`(
  PTRACCAT_AmIndian   = as.integer(pr == "am indian/alaskan"),
  PTRACCAT_Asian      = as.integer(pr == "asian"),
  PTRACCAT_Black      = as.integer(pr == "black"),
  PTRACCAT_MoreOne    = as.integer(pr == "more than one"),
  PTRACCAT_HawaiianPI = as.integer(pr == "hawaiian/other pi")
)]
# (If "white", all dummies above are 0.)

# Marital status (PTMARRY strings; Married is reference)
pm <- norm_lower(m$PTMARRY)
m[, `:=`(
  PTMARRY_Divorced     = as.integer(pm == "divorced"),
  PTMARRY_Widowed      = as.integer(pm == "widowed"),
  PTMARRY_NeverMarried = as.integer(pm == "never married")
)]
# ("married" → all zeros)

# Cohort / Phase (COLPROT strings; ADNI1 is reference)
cp <- toupper(str_trim(m$COLPROT))
m[, `:=`(
  COHORT_ADNIGO = fifelse(cp == "ADNIGO", 1L, fifelse(cp %in% c("ADNI1","ADNI2"), 0L, NA_integer_)),
  COHORT_ADNI2  = fifelse(cp == "ADNI2",  1L, fifelse(cp %in% c("ADNI1","ADNIGO"), 0L, NA_integer_))
)]
# (ADNI1 → both 0; anything else → NA)

# --- Check duplicate PTIDs across cohorts (just in case) ---
dup_ids <- m[, .N, by = .(PTID_like)][N > 1]
if (nrow(dup_ids) > 0) {
  cat("\n== Duplicate PTIDs across cohorts ==\n")
  print(dup_ids)
} else {
  cat("\nNo duplicate PTIDs across cohorts.\n")
}

# --- Distribution of cohorts (COLPROT) ---
cat("\n== Cohort distribution (COLPROT) ==\n")
print(m[, .N, by = COLPROT][order(-N)])

# ========= 6) Merge PCs (from plink2 --pca 10) =========
pcs <- fread(evec_path, header = FALSE)
setnames(pcs, c("FID","IID", paste0("PC", seq_len(ncol(pcs)-2))))
pc_keep <- paste0("PC", 1:10)
pcs <- pcs[, c("FID","IID", pc_keep), with = FALSE]
m <- merge(m, pcs, by = c("FID","IID"), all.x = TRUE)

# ========= 7) Write files for PLINK =========
# Phenotype (from m to keep alignment)
fwrite(m[, .(FID, IID, DIAG01)],
       "ADNI_pheno_ADvsCN1.txt", sep = "\t", quote = FALSE, na = "NA")

# Covariates
covar_cols <- c(
  "FID","IID",
  "AGE","PTEDUCAT","PTGENDER_Female",
  "PTRACCAT_AmIndian","PTRACCAT_Asian","PTRACCAT_Black","PTRACCAT_MoreOne","PTRACCAT_HawaiianPI",
  "PTETHCAT_HispLatino",
  "PTMARRY_Divorced","PTMARRY_Widowed","PTMARRY_NeverMarried",
  "COHORT_ADNIGO","COHORT_ADNI2",
  pc_keep
)
covar_cols <- intersect(covar_cols, names(m))
covar <- m[, ..covar_cols]
fwrite(covar, "ADNI_covar_ADvsCN1.txt", sep = "\t", quote = FALSE, na = "NA")

# Keep list: DIAG01 ∈ {1,2} AND present in covar
valid_ids <- covar[, .(FID, IID)]
keep <- merge(m[, .(FID, IID, DIAG01)], valid_ids, by = c("FID","IID"))
keep <- keep[DIAG01 %in% c(1L, 2L), .(FID, IID)]
fwrite(keep, "ADNI_keep_ADvsCN1.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# ========= 8) Quick summaries =========
cat("\n== Final DIAG01 distribution (from m) ==\n")
print(m[, .N, by = DIAG01][order(DIAG01)])

cat("\n== Covariate missingness (excluding FID/IID) ==\n")
print(covar[, lapply(.SD, function(x) sum(is.na(x))),
            .SDcols = setdiff(names(covar), c("FID","IID"))])

cat("\n== Cohort dummy counts (ADNI1 is reference → both zeros) ==\n")
for (col in c("COHORT_ADNIGO","COHORT_ADNI2")) {
  cat("\n--", col, "--\n")
  print(covar[, .N, by = get(col)][order(-N)])
}

# --- Rows in merged table with ANY missing pheno/covar (useful debug) ---
cols_check <- c("DIAG01","AGE","PTEDUCAT","PTGENDER_Female",
                "PTRACCAT_AmIndian","PTRACCAT_Asian","PTRACCAT_Black",
                "PTRACCAT_MoreOne","PTRACCAT_HawaiianPI",
                "PTETHCAT_HispLatino",
                "PTMARRY_Divorced","PTMARRY_Widowed","PTMARRY_NeverMarried",
                "COHORT_ADNIGO","COHORT_ADNI2",
                paste0("PC", 1:10))
cols_check <- intersect(cols_check, names(m))
m[, n_missing := rowSums(is.na(.SD)), .SDcols = cols_check]
rows_missing <- m[n_missing > 0, .(FID, IID, n_missing)]
cat("\n== Rows with ANY missing phenotype/covariate ==\n")
print(rows_missing)
cat("\nTotal with missingness: ", nrow(rows_missing),
    " / ", nrow(m), " (",
    sprintf('%.1f', 100*nrow(rows_missing)/nrow(m)), "%)\n", sep = "")

# Save details
fwrite(rows_missing, "adnimerge_rows_with_missing.tsv", sep = "\t")

cat("\nWrote:\n  ADNI_pheno_ADvsCN1.txt\n  ADNI_covar_ADvsCN1.txt\n  ADNI_keep_ADvsCN1.txt\n")