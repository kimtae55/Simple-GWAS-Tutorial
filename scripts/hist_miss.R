# robust_hist_miss.R

read_missing <- function(path, type=c("sample","variant")){
  type <- match.arg(type)
  # try with header=TRUE first
  df <- try(read.table(path, header=TRUE, sep="", check.names=FALSE), silent=TRUE)
  bad_header <- inherits(df, "try-error") ||
    all(grepl("^[0-9._-]+$", names(df))) ||  # names look like numeric row
    any(duplicated(names(df)))               # duplicated "headers"
  if (bad_header) {
    df <- read.table(path, header=FALSE, sep="", check.names=FALSE)
    # heuristic column names (PLINK 2 typical layout = 5 cols)
    if (ncol(df) == 5) {
      colnames(df) <- if (type == "sample")
        c("FID","IID","N_MISS","N_CALLED","F_MISS") else
          c("CHR","ID","N_MISS","N_CALLED","F_MISS")
    } else {
      # generic names V1..Vn; last col = fraction missing
      colnames(df) <- paste0("V", seq_len(ncol(df)))
      colnames(df)[ncol(df)] <- "F_MISS"
    }
  }
  df
}

pick_fraction <- function(df){
  nm <- names(df)
  i <- grep("^F(_|\\.)?MISS$", nm, ignore.case=TRUE)
  if (length(i) == 0) i <- ncol(df)  # fallback: last column
  as.numeric(df[[ i[1] ]])
}

# locate files
sm_path <- if (file.exists("plink.imiss")) "plink.imiss" else "plink2.smiss"
vm_path <- if (file.exists("plink.lmiss")) "plink.lmiss" else "plink2.vmiss"

indmiss <- read_missing(sm_path, "sample")
snpmiss <- read_missing(vm_path, "variant")

x_i <- pick_fraction(indmiss)
x_v <- pick_fraction(snpmiss)

pdf("histimiss.pdf"); hist(x_i, main="Histogram: individual missingness", xlab="F_MISS (individual)"); dev.off()
pdf("histlmiss.pdf"); hist(x_v, main="Histogram: variant missingness",   xlab="F_MISS (variant)");   dev.off()

cat("OK: wrote histimiss.pdf and histlmiss.pdf\n")