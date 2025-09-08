# gender_check.R (robust)

# Read with whitespace delimiter; keep header; don't mangle names
gender <- read.table(
  "plink2.sexcheck",
  header = TRUE,
  sep = "",
  na.strings = c("NA", "."),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""   # don't treat '#' lines as comments
)

# Normalize column names: strip leading '#', trim spaces
nm <- colnames(gender)
nm <- sub("^#", "", nm)
nm <- trimws(nm)
# make duplicates unique but preserve visible names
nm <- make.unique(nm, sep = "_")
colnames(gender) <- nm

# Helper: find the F column robustly
find_F_col <- function(nm) {
  # exact "F"
  i <- which(toupper(nm) == "F")
  if (length(i) == 0) {
    # tolerate variants like "F." "F_MISS" (rare in sexcheck), or stray spaces
    i <- grep("^F(\\b|\\.|_|$)", nm, ignore.case = TRUE)
  }
  if (length(i) == 0) stop("Could not find 'F' column. Names: ", paste(nm, collapse=", "))
  i[1]
}

iF <- find_F_col(colnames(gender))

# Coerce needed columns to numeric/integer safely
gender[[iF]] <- suppressWarnings(as.numeric(gender[[iF]]))

# PEDSEX can be read as char; coerce to integer
if ("PEDSEX" %in% names(gender)) {
  gender$PEDSEX <- suppressWarnings(as.integer(gender$PEDSEX))
}

# Optional: YRATE if present
if ("YRATE" %in% names(gender)) {
  gender$YRATE <- suppressWarnings(as.numeric(gender$YRATE))
}

# Plot PDFs
pdf("Gender_check.pdf")
hist(gender[[iF]], main = "Gender", xlab = "F")
dev.off()

pdf("Men_check.pdf")
male <- subset(gender, !is.na(PEDSEX) & PEDSEX == 1)
hist(male[[iF]], main = "Men", xlab = "F")
dev.off()

pdf("Women_check.pdf")
female <- subset(gender, !is.na(PEDSEX) & PEDSEX == 2)
hist(female[[iF]], main = "Women", xlab = "F")
dev.off()

# Optional: print names to confirm
cat("Columns:", paste(colnames(gender), collapse=", "), "\n")