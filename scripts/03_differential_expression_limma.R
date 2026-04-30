# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 3: Differential expression analysis using limma
# ============================================================

# Purpose:
# This script compares Cancer vs Normal samples using limma.
# It creates a differential expression results table.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(limma)
})

# -----------------------------
# 1. Read expression matrix and sample groups
# -----------------------------
expr_df <- read_csv("data/processed/GSE15852_expression_matrix.csv", show_col_types = FALSE)
sample_groups <- read_csv("data/processed/GSE15852_sample_groups.csv", show_col_types = FALSE)
feature_annot <- read_csv("data/processed/GSE15852_feature_annotation.csv", show_col_types = FALSE)

# Convert expression table to matrix
expr_matrix <- expr_df %>%
  column_to_rownames("probe_id") %>%
  as.matrix()

# -----------------------------
# 2. Match expression columns to sample metadata
# -----------------------------
common_samples <- intersect(colnames(expr_matrix), sample_groups$sample_id)

cat("Samples in expression matrix:", ncol(expr_matrix), "\n")
cat("Samples in group table:", nrow(sample_groups), "\n")
cat("Matched samples:", length(common_samples), "\n")

expr_matrix <- expr_matrix[, common_samples]
sample_groups <- sample_groups %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

# Confirm order is correct
stopifnot(all(colnames(expr_matrix) == sample_groups$sample_id))

# -----------------------------
# 3. Check whether log2 transform is needed
# -----------------------------
expr_range <- range(expr_matrix, na.rm = TRUE)
cat("Expression value range before transformation:", expr_range[1], "to", expr_range[2], "\n")

# Microarray data may be already log2-transformed.
# If values are very large, we apply log2 transformation.
if (expr_range[2] > 100) {
  cat("Large expression values detected. Applying log2(x + 1) transformation.\n")
  expr_matrix <- log2(expr_matrix + 1)
} else {
  cat("Expression values appear already log-scale. No log2 transformation applied.\n")
}

expr_range_after <- range(expr_matrix, na.rm = TRUE)
cat("Expression value range after check:", expr_range_after[1], "to", expr_range_after[2], "\n")

# Save the analysis-ready expression matrix
expr_ready <- as.data.frame(expr_matrix) %>%
  rownames_to_column("probe_id")
write_csv(expr_ready, "data/processed/GSE15852_expression_matrix_analysis_ready.csv")

# -----------------------------
# 4. Create limma design matrix
# -----------------------------
sample_groups$group <- factor(sample_groups$group, levels = c("Normal", "Cancer"))

design <- model.matrix(~ 0 + group, data = sample_groups)
colnames(design) <- levels(sample_groups$group)

cat("\nDesign matrix preview:\n")
print(head(design))

# -----------------------------
# 5. Fit limma model
# -----------------------------
fit <- lmFit(expr_matrix, design)

contrast_matrix <- makeContrasts(
  Cancer_vs_Normal = Cancer - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 6. Extract differential expression results
# -----------------------------
de_results <- topTable(
  fit2,
  coef = "Cancer_vs_Normal",
  number = Inf,
  adjust.method = "BH"
) %>%
  rownames_to_column("probe_id")

# -----------------------------
# 7. Add gene annotation if available
# -----------------------------
# GPL96 often contains columns such as Gene Symbol, Gene title, Entrez Gene.
# We keep all annotation columns to preserve useful information.
de_annotated <- de_results %>%
  left_join(feature_annot, by = "probe_id")

# -----------------------------
# 8. Add significance labels
# -----------------------------
de_annotated <- de_annotated %>%
  mutate(
    regulation = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Upregulated in Cancer",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated in Cancer",
      TRUE ~ "Not significant"
    )
  )

# -----------------------------
# 9. Save results
# -----------------------------
write_csv(de_annotated, "results/differential_expression_results.csv")

summary_counts <- de_annotated %>%
  count(regulation)

write_csv(summary_counts, "results/03_differential_expression_summary.csv")

cat("\nDifferential expression summary:\n")
print(summary_counts)

cat("\nTop 10 differential expression results:\n")
print(head(de_annotated, 10))

cat("\nFiles created:\n")
cat("data/processed/GSE15852_expression_matrix_analysis_ready.csv\n")
cat("results/differential_expression_results.csv\n")
cat("results/03_differential_expression_summary.csv\n")
