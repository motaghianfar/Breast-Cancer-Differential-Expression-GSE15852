# Differential Gene Expression in Breast Cancer vs Normal Tissue

## Project Overview

This beginner-friendly bioinformatics project identifies genes that are differentially expressed between breast cancer tissue and normal breast tissue. The goal is to find genes that are significantly upregulated or downregulated in cancer. This type of analysis is commonly used in transcriptomics, cancer genomics, biomarker discovery, and disease mechanism research.

## Research Topic

**Exact topic:** Identification of differentially expressed genes between breast cancer and normal breast tissue using gene expression data.

**Objective:** Find genes significantly upregulated or downregulated in cancer.

**Bioinformatics area:** Transcriptomics / Cancer Genomics

## Dataset

- **Dataset accession:** GSE15852
- **Repository:** NCBI Gene Expression Omnibus GEO
- **Platform:** GPL96 Affymetrix Human Genome U133A Array
- **Data type:** Gene expression
- **Samples:** 86 total
  - 43 breast cancer samples
  - 43 normal breast tissue samples
- **Features tested:** 22,283 probes

## Method Summary

This project uses the `GEOquery` R package to download the dataset from GEO, then uses `limma` to compare expression between cancer and normal samples. Since GSE15852 provides processed Affymetrix microarray expression data, `limma` is the appropriate method for differential expression analysis. Genes/probes were considered significant if they had:

- Adjusted P-value < 0.05
- Absolute log2 fold change >= 1

## Key Results

| Category | Count |
|---|---:|
| Upregulated in cancer | 141 |
| Downregulated in cancer | 267 |
| Not significant | 21,875 |
| Total probes tested | 22,283 |

## Example Top Upregulated Genes

Examples of genes upregulated in cancer include:

- CD24
- KRT19
- EPCAM
- TACSTD2
- RGS1

## Example Top Downregulated Genes

Examples of genes downregulated in cancer include:

- PPP1R1A
- RBP4
- PDE3B
- GYG2
- ANGPTL4
- MAOA

## Figures

The project generates the following figures:

1. **Volcano plot**  
   `figures/volcano_plot.png`

2. **PCA plot**  
   `figures/pca_plot.png`

3. **Heatmap of top 20 differentially expressed genes**  
   `figures/heatmap_top20_genes.png`

## Results Tables

Important result files are saved in the `results/` folder:

- `results/differential_expression_results.csv`
- `results/project_summary.csv`
- `results/top20_differentially_expressed_genes.csv`
- `results/top10_upregulated_genes.csv`
- `results/top10_downregulated_genes.csv`
- `results/pca_coordinates.csv`

## Project Structure

```text
Project1/
├── data/
│   ├── raw/
│   └── processed/
├── figures/
├── results/
├── scripts/
├── report/
├── docs/
├── README.md
├── .gitignore
└── commands.shcat > README.md <<'EOF'
# Differential Gene Expression in Breast Cancer vs Normal Tissue

## Project Overview

This beginner-friendly bioinformatics project identifies genes that are differentially expressed between breast cancer tissue and normal breast tissue. The goal is to find genes that are significantly upregulated or downregulated in cancer. This type of analysis is commonly used in transcriptomics, cancer genomics, biomarker discovery, and disease mechanism research.

## Research Topic

**Exact topic:** Identification of differentially expressed genes between breast cancer and normal breast tissue using gene expression data.

**Objective:** Find genes significantly upregulated or downregulated in cancer.

**Bioinformatics area:** Transcriptomics / Cancer Genomics

## Dataset

- **Dataset accession:** GSE15852
- **Repository:** NCBI Gene Expression Omnibus GEO
- **Platform:** GPL96 Affymetrix Human Genome U133A Array
- **Data type:** Gene expression
- **Samples:** 86 total
  - 43 breast cancer samples
  - 43 normal breast tissue samples
- **Features tested:** 22,283 probes

## Method Summary

This project uses the `GEOquery` R package to download the dataset from GEO, then uses `limma` to compare expression between cancer and normal samples. Since GSE15852 provides processed Affymetrix microarray expression data, `limma` is the appropriate method for differential expression analysis. Genes/probes were considered significant if they had:

- Adjusted P-value < 0.05
- Absolute log2 fold change >= 1

## Key Results

| Category | Count |
|---|---:|
| Upregulated in cancer | 141 |
| Downregulated in cancer | 267 |
| Not significant | 21,875 |
| Total probes tested | 22,283 |

## Example Top Upregulated Genes

Examples of genes upregulated in cancer include:

- CD24
- KRT19
- EPCAM
- TACSTD2
- RGS1

## Example Top Downregulated Genes

Examples of genes downregulated in cancer include:

- PPP1R1A
- RBP4
- PDE3B
- GYG2
- ANGPTL4
- MAOA

## Figures

The project generates the following figures:

1. **Volcano plot**  
   `figures/volcano_plot.png`

2. **PCA plot**  
   `figures/pca_plot.png`

3. **Heatmap of top 20 differentially expressed genes**  
   `figures/heatmap_top20_genes.png`

## Results Tables

Important result files are saved in the `results/` folder:

- `results/differential_expression_results.csv`
- `results/project_summary.csv`
- `results/top20_differentially_expressed_genes.csv`
- `results/top10_upregulated_genes.csv`
- `results/top10_downregulated_genes.csv`
- `results/pca_coordinates.csv`

## Project Structure

```text
Project1/
├── data/
│   ├── raw/
│   └── processed/
├── figures/
├── results/
├── scripts/
├── report/
├── docs/
├── README.md
├── .gitignore
└── commands.sh
