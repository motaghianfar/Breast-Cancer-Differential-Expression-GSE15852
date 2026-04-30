# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 4: Volcano plot
# ============================================================

# Purpose:
# This script creates a volcano plot from the limma differential
# expression results. The plot highlights significantly upregulated
# and downregulated genes in breast cancer.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# -----------------------------
# 1. Read differential expression results
# -----------------------------
de_results <- read_csv("results/differential_expression_results.csv", show_col_types = FALSE)

# -----------------------------
# 2. Prepare labels
# -----------------------------
# Use Gene Symbol where available; otherwise use probe_id.
de_results <- de_results %>%
  mutate(
    gene_label = ifelse(
      !is.na(`Gene Symbol`) & `Gene Symbol` != "",
      `Gene Symbol`,
      probe_id
    ),
    neg_log10_adj_p = -log10(adj.P.Val)
  )

# Select top genes for labels:
# strongest significant genes by adjusted p-value
top_labels <- de_results %>%
  filter(regulation != "Not significant") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 15)

# -----------------------------
# 3. Create volcano plot
# -----------------------------
volcano_plot <- ggplot(de_results, aes(x = logFC, y = neg_log10_adj_p)) +
  geom_point(aes(color = regulation), alpha = 0.7, size = 1.6) +
  scale_color_manual(
    values = c(
      "Upregulated in Cancer" = "red",
      "Downregulated in Cancer" = "blue",
      "Not significant" = "gray70"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_label),
    size = 3,
    max.overlaps = 50
  ) +
  labs(
    title = "Volcano Plot: Breast Cancer vs Normal Tissue",
    subtitle = "GSE15852 differential expression analysis using limma",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Regulation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# -----------------------------
# 4. Save figure
# -----------------------------
ggsave(
  filename = "figures/volcano_plot.png",
  plot = volcano_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "figures/volcano_plot.pdf",
  plot = volcano_plot,
  width = 9,
  height = 6
)

cat("Volcano plot created:\n")
cat("figures/volcano_plot.png\n")
cat("figures/volcano_plot.pdf\n")
