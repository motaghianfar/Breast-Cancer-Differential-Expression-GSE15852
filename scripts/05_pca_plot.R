# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 5: PCA plot
# ============================================================

# Purpose:
# This script creates a PCA plot to show whether Cancer and Normal
# samples separate based on overall gene expression patterns.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# -----------------------------
# 1. Read analysis-ready expression matrix and sample groups
# -----------------------------
expr_df <- read_csv("data/processed/GSE15852_expression_matrix_analysis_ready.csv", show_col_types = FALSE)
sample_groups <- read_csv("data/processed/GSE15852_sample_groups.csv", show_col_types = FALSE)

expr_matrix <- expr_df %>%
  column_to_rownames("probe_id") %>%
  as.matrix()

# -----------------------------
# 2. Match sample order
# -----------------------------
common_samples <- intersect(colnames(expr_matrix), sample_groups$sample_id)

expr_matrix <- expr_matrix[, common_samples]
sample_groups <- sample_groups %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

stopifnot(all(colnames(expr_matrix) == sample_groups$sample_id))

# -----------------------------
# 3. Select most variable probes
# -----------------------------
# PCA is clearer when we use probes that vary across samples.
# Here we use the top 1000 most variable probes.
probe_variance <- apply(expr_matrix, 1, var, na.rm = TRUE)

top_variable_probes <- names(sort(probe_variance, decreasing = TRUE))[1:1000]

expr_top <- expr_matrix[top_variable_probes, ]

# -----------------------------
# 4. Run PCA
# -----------------------------
# prcomp expects samples as rows and genes/probes as columns,
# so we transpose the expression matrix.
pca_result <- prcomp(t(expr_top), scale. = TRUE)

# Calculate percent variance explained by PC1 and PC2
pca_var <- pca_result$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

pca_df <- as.data.frame(pca_result$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(sample_groups, by = "sample_id")

# -----------------------------
# 5. Create PCA plot
# -----------------------------
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA Plot: Breast Cancer vs Normal Tissue",
    subtitle = "GSE15852 expression profiles using top 1000 variable probes",
    x = paste0("PC1 (", pca_var_percent[1], "% variance)"),
    y = paste0("PC2 (", pca_var_percent[2], "% variance)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# -----------------------------
# 6. Save PCA figure and PCA coordinates
# -----------------------------
ggsave(
  filename = "figures/pca_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "figures/pca_plot.pdf",
  plot = pca_plot,
  width = 8,
  height = 6
)

write_csv(pca_df, "results/pca_coordinates.csv")

cat("PCA plot created:\n")
cat("figures/pca_plot.png\n")
cat("figures/pca_plot.pdf\n")
cat("results/pca_coordinates.csv\n")

cat("\nVariance explained:\n")
cat("PC1:", pca_var_percent[1], "%\n")
cat("PC2:", pca_var_percent[2], "%\n")
