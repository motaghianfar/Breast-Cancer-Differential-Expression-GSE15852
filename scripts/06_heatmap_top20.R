# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 6: Heatmap of top 20 differentially expressed genes
# ============================================================

# Purpose:
# This script creates a heatmap of the top 20 most significant
# differentially expressed probes/genes between Cancer and Normal samples.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
})

# -----------------------------
# 1. Read required files
# -----------------------------
expr_df <- read_csv("data/processed/GSE15852_expression_matrix_analysis_ready.csv", show_col_types = FALSE)
sample_groups <- read_csv("data/processed/GSE15852_sample_groups.csv", show_col_types = FALSE)
de_results <- read_csv("results/differential_expression_results.csv", show_col_types = FALSE)

expr_matrix <- expr_df %>%
  column_to_rownames("probe_id") %>%
  as.matrix()

# -----------------------------
# 2. Select top 20 significant genes/probes
# -----------------------------
top20 <- de_results %>%
  filter(regulation != "Not significant") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

top20_probes <- top20$probe_id

# -----------------------------
# 3. Extract expression values for top 20
# -----------------------------
heatmap_matrix <- expr_matrix[top20_probes, ]

# -----------------------------
# 4. Use gene symbols as row labels where possible
# -----------------------------
gene_labels <- ifelse(
  !is.na(top20$`Gene Symbol`) & top20$`Gene Symbol` != "",
  top20$`Gene Symbol`,
  top20$probe_id
)

# Make labels unique in case repeated gene symbols exist
gene_labels <- make.unique(gene_labels)
rownames(heatmap_matrix) <- gene_labels

# -----------------------------
# 5. Match sample annotation
# -----------------------------
sample_groups <- sample_groups %>%
  filter(sample_id %in% colnames(heatmap_matrix)) %>%
  arrange(match(sample_id, colnames(heatmap_matrix)))

stopifnot(all(sample_groups$sample_id == colnames(heatmap_matrix)))

annotation_col <- sample_groups %>%
  select(sample_id, group) %>%
  column_to_rownames("sample_id")

# -----------------------------
# 6. Save top 20 gene table
# -----------------------------
top20_table <- top20 %>%
  select(
    probe_id,
    `Gene Symbol`,
    `Gene Title`,
    logFC,
    P.Value,
    adj.P.Val,
    regulation
  )

write_csv(top20_table, "results/top20_differentially_expressed_genes.csv")

# -----------------------------
# 7. Create heatmap
# -----------------------------
# scale = "row" means each gene is standardized across samples.
# This makes relative high/low expression easier to compare visually.

png("figures/heatmap_top20_genes.png", width = 2400, height = 1800, res = 300)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
  main = "Top 20 Differentially Expressed Genes: Cancer vs Normal"
)

dev.off()

pdf("figures/heatmap_top20_genes.pdf", width = 10, height = 8)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
  main = "Top 20 Differentially Expressed Genes: Cancer vs Normal"
)

dev.off()

cat("Heatmap created:\n")
cat("figures/heatmap_top20_genes.png\n")
cat("figures/heatmap_top20_genes.pdf\n")
cat("results/top20_differentially_expressed_genes.csv\n")

cat("\nTop 20 genes/probes used in heatmap:\n")
print(top20_table)
