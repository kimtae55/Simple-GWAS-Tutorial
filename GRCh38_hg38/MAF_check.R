# MAF_check.R  (robust for PLINK 2 .afreq, mono- and multi-allelic)
af <- read.table(
  "MAF_check.afreq",
  header = TRUE,
  sep = "\t",
  comment.char = "",         # keep '#CHROM' header
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Ensure ALT_FREQS exists and is character
if (!"ALT_FREQS" %in% names(af)) stop("ALT_FREQS column not found.")
af$ALT_FREQS <- as.character(af$ALT_FREQS)

# Parse ALT_FREQS (may be '0', '0.187', or '0.2,0.05' for multi-allelic)
split_freqs <- strsplit(af$ALT_FREQS, ",", fixed = TRUE)
alt_sums    <- vapply(split_freqs, function(v) sum(as.numeric(v), na.rm = TRUE), numeric(1))
# ref frequency = 1 - sum(alt freqs); clamp to [0,1] to be safe
ref_freq    <- pmax(0, pmin(1, 1 - alt_sums))

# Per row, minor allele frequency = min(ref_freq, each alt freq)
row_maf <- vapply(seq_along(split_freqs), function(i) {
  alts <- suppressWarnings(as.numeric(split_freqs[[i]]))
  if (length(alts) == 0 || all(is.na(alts))) alts <- 0  # handles ALT="."
  min(c(ref_freq[i], alts), na.rm = TRUE)
}, numeric(1))

# Histogram
pdf("MAF_distribution.pdf")
hist(row_maf, breaks = 50, main = "MAF distribution", xlab = "Minor Allele Frequency")
dev.off()

cat("OK. Variants:", nrow(af), "  Wrote: MAF_distribution.pdf\n")