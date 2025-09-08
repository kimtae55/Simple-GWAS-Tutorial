# ================================
# Load PCA results from merged data
# ================================
# File from PLINK2: merged_ADNI_1KG_PCA.eigenvec
ev <- read.table("merged_ADNI_1KG_PCA.eigenvec", header = FALSE)
colnames(ev) <- c("FID", "IID", paste0("PC", 1:(ncol(ev) - 2)))

# ================================
# Load 1000 Genomes population labels
# ================================
# Download this first if you haven't:
# wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20131219.populations.tsv
pop <- read.table("20131219.populations.tsv", header = TRUE, sep = "\t")

# Merge PCA data with 1000G super-population labels
ref <- merge(ev, pop, by.x = "IID", by.y = "Sample", all.x = TRUE)

# Add a column to distinguish ADNI vs 1000G reference samples
ref$Dataset <- ifelse(is.na(ref$Super_Population), "ADNI", "1KG")

# ================================
# Set color palette for super-populations
# ================================
cols <- c(
  AFR = "#E41A1C",  # African - red
  EUR = "#377EB8",  # European - blue
  EAS = "#4DAF4A",  # East Asian - green
  AMR = "#FF7F00",  # Admixed American - orange
  SAS = "#984EA3"   # South Asian - purple
)

# Assign gray to ADNI samples for initial plotting
ref$Color <- ifelse(ref$Dataset == "ADNI", "black", cols[ref$Super_Population])

# ================================
# Plot PCA: 1000G populations + ADNI overlay
# ================================
png("ADNI_vs_1KG_PCA.png", width = 1200, height = 900, res = 120)

plot(ref$PC1, ref$PC2,
     col = ref$Color,
     pch = ifelse(ref$Dataset == "ADNI", 17, 19), # triangles = ADNI
     cex = ifelse(ref$Dataset == "ADNI", 1.0, 0.7),
     xlab = "PC1", ylab = "PC2",
     main = "ADNI Samples Projected onto 1000 Genomes")

# Add legend for super-populations + ADNI
legend("topright",
       legend = c(names(cols), "ADNI"),
       col = c(cols, "black"),
       pch = c(rep(19, length(cols)), 17),
       bty = "n")

grid()
dev.off()

# ================================
# Save top 10 PCs for ADNI samples
# ================================
write.csv(ref[ref$Dataset == "ADNI", c("FID", "IID", paste0("PC", 1:10))],
          "ADNI_top10_PCs.csv",
          row.names = FALSE)