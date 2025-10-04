# hwe.R  (PLINK 2 --hardy output)

hwe <- read.table(
  "plink2.hardy",
  header = TRUE,
  sep = "",              # any whitespace
  comment.char = "",     # KEEP '#CHROM' header
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# fix header name '#CHROM' → 'CHROM'
colnames(hwe) <- sub("^#", "", colnames(hwe))

# coerce P to numeric (it may be read as char)
hwe$P <- suppressWarnings(as.numeric(hwe$P))

# sanity check (optional)
# print(colnames(hwe)); str(hwe$P)

pdf("histhwe.pdf")
hist(hwe$P, main="Histogram of HWE P-values", xlab="HWE P-value", breaks=50)
dev.off()

hwe_zoom <- subset(hwe, P < 1e-5)
pdf("histhwe_below_threshold.pdf")
hist(hwe_zoom$P, main="Histogram of SNPs deviating from HWE (P < 1e-5)", xlab="HWE P-value", breaks=50)
dev.off()

write.table(hwe_zoom, "plinkzoomhwe.hardy", sep="\t", quote=FALSE, row.names=FALSE)