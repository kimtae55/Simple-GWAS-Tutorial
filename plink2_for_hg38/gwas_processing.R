suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ========= Paths (edit if needed) =========
fam_path   <- "/Users/taehyo/Local_Document/data_l0ipls/ADNI_MERGED_FINAL/ADNI_qc_final.fam"
dxsum_path <- "/Users/taehyo/Local_Document/data_l0ipls/clinical/DXSUM_11Aug2025.csv"
ptdem_path <- "/Users/taehyo/Local_Document/data_l0ipls/clinical/PTDEMOG_28Jul2025.csv"
evec_path  <- "/Users/taehyo/Local_Document/data_l0ipls/ADNI_MERGED_FINAL/ADNI_pca.eigenvec"  # from: plink2 --pca 10

# ========= Helpers =========
parse_date_flex <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(as.IDate(NA), length(x))
  # YYYY-MM-DD or YYYY/MM/DD
  iso <- grepl("^\\d{4}[-/]\\d{2}[-/]\\d{2}$", x)
  out[iso] <- as.IDate(gsub("/", "-", x[iso]))
  # MDY / DMY with / or -
  mdy_dmy <- is.na(out) & grepl("^\\d{1,2}[/-]\\d{1,2}[/-]\\d{4}$", x)
  if (any(mdy_dmy)) {
    fmts <- c("%m/%d/%Y","%m-%d-%Y","%d/%m/%Y","%d-%m-%Y")
    tmp <- rep(NA_real_, sum(mdy_dmy))
    for (fmt in fmts) {
      hit <- is.na(tmp)
      if (any(hit)) {
        parsed <- as.Date(x[mdy_dmy][hit], format = fmt)
        tmp[hit] <- as.numeric(parsed)
      }
    }
    out[mdy_dmy] <- as.IDate(as.Date(tmp, origin = "1970-01-01"))
  }
  # MM/YYYY → day 15
  mY <- is.na(out) & grepl("^\\d{1,2}[/-]\\d{4}$", x)
  if (any(mY)) {
    m <- as.integer(sub("[/\\-].*$", "", x[mY])); y <- as.integer(sub("^.*[/\\-]", "", x[mY]))
    m[is.na(m) | m < 1 | m > 12] <- 6L
    out[mY] <- as.IDate(sprintf("%04d-%02d-15", y, m))
  }
  # YYYY → July 1
  Y <- is.na(out) & grepl("^\\d{4}$", x)
  if (any(Y)) {
    y <- as.integer(x[Y])
    out[Y] <- as.IDate(sprintf("%04d-07-01", y))
  }
  out
}
age_years <- function(dob, refdate) as.numeric((refdate - dob) / 365.25)
to_chr <- function(x) { x <- as.character(x); x[is.na(x)] <- NA_character_; trimws(x) }
has_code <- function(x, code) {
  x <- to_chr(x)
  ifelse(is.na(x) | x == "", FALSE,
         grepl(paste0("(^|\\|)", code, "($|\\|)"), x))
}

# ========= 1) FAM (SEX from fam; PTID_like) =========
fam <- fread(fam_path, header = FALSE,
             col.names = c("FID","IID","MID","PID","SEX","PHENO"))
fam[, PTID_like := toupper(sub("^[^_]+_", "", IID))]   # e.g., "1_014_S_0520" -> "014_S_0520"
# Female dummy from fam: 1=male, 2=female
fam[, PTGENDER_Female := fifelse(SEX == 2, 1L, fifelse(SEX == 1, 0L, NA_integer_))]

# ========= 2) DXSUM baseline (VISCODE2 == 'bl'), DIAGNOSIS -> DIAG01 =========
dx <- fread(dxsum_path)
stopifnot(all(c("PTID","VISCODE2","DIAGNOSIS","EXAMDATE") %in% names(dx)))
dx[, PTID := toupper(str_trim(PTID))]
dx_bl <- dx[tolower(VISCODE2) == "bl"]
dx_bl[, DIAGNOSIS := suppressWarnings(as.integer(DIAGNOSIS))]
dx_bl[, DIAG01 := fifelse(DIAGNOSIS == 1L, 1L,    # CN
                          fifelse(DIAGNOSIS == 3L, 2L,   # AD
                                  -9L))]                 # others
dx_bl[, EXAMDATE := parse_date_flex(EXAMDATE)]
dx_bl <- dx_bl[, .(PTID, DIAG01, EXAMDATE)]
dx_bl <- dx_bl[!duplicated(PTID)]  # keep one baseline row per PTID

# Anchor to IDs present in both FAM and DXSUM baseline
anchor <- merge(
  fam[, .(FID, IID, PTID_like, PTGENDER_Female)],
  dx_bl[, .(PTID, DIAG01, EXAMDATE)],
  by.x = "PTID_like",   # length 1
  by.y = "PTID",        # length 1
  all = FALSE
)

# ========= 3) PTDEMOG (for AGE, PTEDUCAT, race, ethnicity, marital) =========
pd <- fread(ptdem_path)
stopifnot(all(c("PTID","PTDOB","PTEDUCAT","PTETHCAT","PTRACCAT","PTMARRY") %in% names(pd)))
pd[, PTID := toupper(str_trim(PTID))]
pd[, PTDOB_DATE := parse_date_flex(PTDOB)]
pd <- pd[, .(PTID, PTDOB_DATE, PTEDUCAT, PTETHCAT, PTRACCAT, PTMARRY)]
pd <- pd[!duplicated(PTID)]

# Merge & compute AGE
m <- merge(anchor, pd, by.x = "PTID_like", by.y = "PTID", all.x = TRUE)
m[, AGE := age_years(PTDOB_DATE, EXAMDATE)]
m[, PTEDUCAT := suppressWarnings(as.numeric(PTEDUCAT))]

# ========= 4) PCs from plink2 --pca 10 =========
pcs <- fread(evec_path, header = FALSE)
setnames(pcs, c("FID","IID", paste0("PC", seq_len(ncol(pcs)-2))))
pc_keep <- paste0("PC", 1:10)
pcs <- pcs[, c("FID","IID", pc_keep), with = FALSE]
m <- merge(m, pcs, by = c("FID","IID"), all.x = TRUE)

# ========= 5) Dummies per your codebooks =========
# Ethnicity (PTETHCAT): 2 = Not Hisp/Latino, 1 = Hisp/Latino
m[, PTETHCAT := to_chr(PTETHCAT)]
m[, PTETHCAT_HispLatino := fifelse(PTETHCAT == "1", 1L,
                                   fifelse(PTETHCAT == "2", 0L, NA_integer_))]

# Race (PTRACCAT) — multi-select strings; White=5 is reference (no dummy)
# Codes: 1=Am Indian/Alaskan, 2=Asian, 3=Hawaiian/Other PI, 4=Black, 5=White, 6=More than one
m[, PTRACCAT := to_chr(PTRACCAT)]
m[, `:=`(
  PTRACCAT_AmIndian   = as.integer(has_code(PTRACCAT, "1")),
  PTRACCAT_Asian      = as.integer(has_code(PTRACCAT, "2")),
  PTRACCAT_HawaiianPI = as.integer(has_code(PTRACCAT, "3")),
  PTRACCAT_Black      = as.integer(has_code(PTRACCAT, "4")),
  PTRACCAT_MoreOne    = as.integer(has_code(PTRACCAT, "6"))
)]
# (White=5 => all the above are 0)

# Marital status (PTMARRY) — Married=1 is reference (no dummy for married)
# Coding: 1=Married, 3=Divorced, 2=Widowed, 4=Never married
m[, PTMARRY := to_chr(PTMARRY)]
m[, `:=`(
  PTMARRY_Divorced     = fifelse(PTMARRY == "3", 1L, 0L),
  PTMARRY_Widowed      = fifelse(PTMARRY == "2", 1L, 0L),
  PTMARRY_NeverMarried = fifelse(PTMARRY == "4", 1L, 0L)
)]

# ========= 6) Write phenotype, covariates, keep list =========
# Phenotype file
fwrite(m[, .(FID, IID, DIAG01)],
       "ADNI_pheno_ADvsCN.txt", sep = "\t", quote = FALSE, na = "NA")

# Covariates (add the PTMARRY dummies; married is reference)
covar_cols <- c(
  "FID","IID",
  "AGE","PTEDUCAT","PTGENDER_Female",
  "PTRACCAT_AmIndian","PTRACCAT_Asian","PTRACCAT_Black","PTRACCAT_MoreOne","PTRACCAT_HawaiianPI",
  "PTETHCAT_HispLatino",
  "PTMARRY_Divorced","PTMARRY_Widowed","PTMARRY_NeverMarried",
  paste0("PC", 1:10)
)
covar_cols <- intersect(covar_cols, names(m))
covar <- m[, ..covar_cols]
fwrite(covar, "ADNI_covar_ADvsCN.txt", sep = "\t", quote = FALSE, na = "NA")

# Keep list: DIAG01 in {1,2} AND present in covar table
valid_ids <- covar[, .(FID, IID)]
keep <- merge(m[, .(FID, IID, DIAG01)], valid_ids, by = c("FID","IID"))
keep <- keep[DIAG01 %in% c(1L, 2L), .(FID, IID)]
fwrite(keep, "ADNI_keep_ADvsCN.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# ========= 7) Quick sanity summaries =========
cat("\n== DIAG01 distribution (all baseline) ==\n")
print(m[, .N, by = DIAG01][order(DIAG01)])

cat("\n== Covariate missingness (excluding FID/IID) ==\n")
print(covar[, lapply(.SD, function(x) sum(is.na(x))),
            .SDcols = setdiff(names(covar), c("FID","IID"))])

cat("\n== Selected dummy counts (1 vs 0) ==\n")
for (col in c("PTGENDER_Female","PTETHCAT_HispLatino",
              "PTRACCAT_AmIndian","PTRACCAT_Asian","PTRACCAT_Black","PTRACCAT_MoreOne","PTRACCAT_HawaiianPI",
              "PTMARRY_Divorced","PTMARRY_Widowed","PTMARRY_NeverMarried")) {
  cat("\n--", col, "--\n")
  print(covar[, .N, by = get(col)][order(-N)])
}

cat("\nWrote:\n  ADNI_pheno_ADvsCN.txt\n  ADNI_covar_ADvsCN.txt\n  ADNI_keep_ADvsCN.txt\n")