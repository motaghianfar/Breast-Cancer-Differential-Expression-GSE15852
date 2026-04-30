# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 7: Create clean summary tables
# ============================================================

# Purpose:
# This script creates beginner-friendly summary tables for the
# README, final report, and GitHub repository.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# -----------------------------
# 1. Read analysis outputs
# -----------------------------
de_results <- read_csv("results/differential_expression_results.csv", show_col_types = FALSE)
sample_groups <- read_csv("data/processed/GSE15852_sample_groups.csv", show_col_types = FALSE)
top20 <- read_csv("results/top20_differentially_expressed_genes.csv", show_col_types = FALSE)

# -----------------------------
# 2. Count samples
# -----------------------------
sample_summary <- sample_groups %>%
  count(group, name = "sample_count")

write_csv(sample_summary, "results/sample_summary.csv")

# -----------------------------
# 3. Count differential expression categories
# -----------------------------
de_summary <- de_results %>%
  count(regulation, name = "probe_count")

write_csv(de_summary, "results/differential_expression_summary_clean.csv")

# -----------------------------
# 4. Create one-row project summary
# -----------------------------
project_summary <- tibble::tibble(
  dataset = "GSE15852",
  platform = "GPL96 Affymetrix Human Genome U133A Array",
  total_samples = nrow(sample_groups),
  normal_samples = sum(sample_groups$group == "Normal"),
  cancer_samples = sum(sample_groups$group == "Cancer"),
  total_probes_tested = nrow(de_results),
  upregulated_in_cancer = sum(de_results$regulation == "Upregulated in Cancer"),
  downregulated_in_cancer = sum(de_results$regulation == "Downregulated in Cancer"),
  significance_cutoff = "Adjusted P-value < 0.05 and |log2FC| >= 1",
  method = "limma"
)

write_csv(project_summary, "results/project_summary.csv")

# -----------------------------
# 5. Create top 10 upregulated and downregulated tables
# -----------------------------
top10_up <- de_results %>%
  filter(regulation == "Upregulated in Cancer") %>%
  arrange(adj.P.Val) %>%
  select(probe_id, `Gene Symbol`, `Gene Title`, logFC, P.Value, adj.P.Val, regulation) %>%
  slice_head(n = 10)

top10_down <- de_results %>%
  filter(regulation == "Downregulated in Cancer") %>%
  arrange(adj.P.Val) %>%
  select(probe_id, `Gene Symbol`, `Gene Title`, logFC, P.Value, adj.P.Val, regulation) %>%
  slice_head(n = 10)

write_csv(top10_up, "results/top10_upregulated_genes.csv")
write_csv(top10_down, "results/top10_downregulated_genes.csv")

# -----------------------------
# 6. Print summary
# -----------------------------
cat("Project summary:\n")
print(project_summary)

cat("\nSample summary:\n")
print(sample_summary)

cat("\nDifferential expression summary:\n")
print(de_summary)

cat("\nFiles created:\n")
cat("results/project_summary.csv\n")
cat("results/sample_summary.csv\n")
cat("results/differential_expression_summary_clean.csv\n")
cat("results/top10_upregulated_genes.csv\n")
cat("results/top10_downregulated_genes.csv\n")
