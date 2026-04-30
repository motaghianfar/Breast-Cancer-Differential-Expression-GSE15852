#!/bin/bash

# ============================================================
# Reproducible commands for Project 1
# Differential Gene Expression in Breast Cancer vs Normal Tissue
# Dataset: GSE15852
# ============================================================

set -e

echo "Running Project 1 analysis workflow..."

Rscript scripts/01_download_inspect_GSE15852.R
Rscript scripts/02_create_sample_groups.R
Rscript scripts/03_differential_expression_limma.R
Rscript scripts/04_volcano_plot.R
Rscript scripts/05_pca_plot.R
Rscript scripts/06_heatmap_top20.R
Rscript scripts/07_create_summary_tables.R

echo "Analysis complete."
echo "Figures are in the figures/ folder."
echo "Tables are in the results/ folder."
