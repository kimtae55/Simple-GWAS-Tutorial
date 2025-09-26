# check_heterozygosity_rate.R  (PLINK 2)
het <- read.table(
  "R_check.het",
  header = TRUE,
  sep = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Clean names: O(HOM)->O_HOM, E(HOM)->E_HOM, drop leading '#'
nm <- colnames(het)
nm <- sub("^#", "", nm)
nm <- gsub("[()]", "_", nm)
colnames(het) <- nm
# Expect: FID, IID, O_HOM, E_HOM, OBS_CT, F

# Coerce numeric
het$O_HOM  <- as.numeric(het$O_HOM)
het$OBS_CT <- as.numeric(het$OBS_CT)
het$F      <- as.numeric(het$F)

# HET_RATE = 1 - O_HOM / OBS_CT
het$HET_RATE <- 1 - (het$O_HOM / het$OBS_CT)

# Plot distribution
pdf("heterozygosity.pdf")
hist(het$HET_RATE, xlab="Heterozygosity Rate", main="Heterozygosity Rate", breaks=50)
dev.off()

# Flag outliers (±3 SD)
mu <- mean(het$HET_RATE, na.rm=TRUE)
sdv <- sd(het$HET_RATE,  na.rm=TRUE)
het_fail <- subset(het, HET_RATE < mu - 3*sdv | HET_RATE > mu + 3*sdv)

# Save full table (for audit)
write.table(het_fail, "fail-het-qc.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Save FID IID only (no header) for PLINK --remove
write.table(het_fail[, c("FID","IID")],
            "het_fail_ind.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)