# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 1: Download and inspect GEO dataset GSE15852
# ============================================================

# Purpose:
# This script downloads the GSE15852 dataset from NCBI GEO,
# extracts the expression matrix and sample metadata,
# and saves them for later analysis.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(readr)
  library(tibble)
  library(dplyr)
})

# -----------------------------
# 1. Create output folders
# -----------------------------
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Download GEO dataset
# -----------------------------
cat("Downloading GSE15852 from GEO...\n")

gse <- getGEO("GSE15852", GSEMatrix = TRUE)

cat("Number of expression sets found:", length(gse), "\n")

# If GEO returns multiple platforms, use the first one for now.
# Later we will inspect whether this is the correct expression set.
eset <- gse[[1]]

# -----------------------------
# 3. Extract expression matrix
# -----------------------------
expr_matrix <- exprs(eset)

cat("Expression matrix dimensions:\n")
cat("Genes/probes:", nrow(expr_matrix), "\n")
cat("Samples:", ncol(expr_matrix), "\n")

# Save expression matrix
expr_df <- as.data.frame(expr_matrix) %>%
  rownames_to_column(var = "probe_id")

write_csv(expr_df, "data/processed/GSE15852_expression_matrix.csv")

# -----------------------------
# 4. Extract sample metadata
# -----------------------------
metadata <- pData(eset) %>%
  rownames_to_column(var = "sample_id")

write_csv(metadata, "data/processed/GSE15852_sample_metadata.csv")

# -----------------------------
# 5. Save feature/gene annotation
# -----------------------------
feature_data <- fData(eset) %>%
  rownames_to_column(var = "probe_id")

write_csv(feature_data, "data/processed/GSE15852_feature_annotation.csv")

# -----------------------------
# 6. Print useful metadata columns
# -----------------------------
cat("\nMetadata columns available:\n")
print(colnames(metadata))

cat("\nFirst few metadata rows:\n")
print(head(metadata[, 1:min(8, ncol(metadata))]))

cat("\nFeature annotation columns available:\n")
print(colnames(feature_data))

cat("\nFirst few feature annotation rows:\n")
print(head(feature_data[, 1:min(8, ncol(feature_data))]))

# -----------------------------
# 7. Save summary
# -----------------------------
summary_text <- c(
  "GSE15852 download and inspection summary",
  paste("Number of expression sets:", length(gse)),
  paste("Genes/probes:", nrow(expr_matrix)),
  paste("Samples:", ncol(expr_matrix)),
  "",
  "Files created:",
  "data/processed/GSE15852_expression_matrix.csv",
  "data/processed/GSE15852_sample_metadata.csv",
  "data/processed/GSE15852_feature_annotation.csv"
)

writeLines(summary_text, "results/01_download_inspect_summary.txt")

cat("\nDone. Files saved in data/processed and results folders.\n")
