# ============================================================
# Project 1: Differential Gene Expression in Breast Cancer
# Step 2: Create tumor vs normal sample group labels
# ============================================================

# Purpose:
# This script reads the GEO metadata and creates a clean sample table.
# Each sample is labeled as either Cancer or Normal based on metadata.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# -----------------------------
# 1. Read sample metadata
# -----------------------------
metadata <- read_csv("data/processed/GSE15852_sample_metadata.csv", show_col_types = FALSE)

# -----------------------------
# 2. Create group labels
# -----------------------------
# We use source_name_ch1 and title because the metadata clearly contains:
# "normal breast tissue" for normal samples
# "breast tumor tissue" for cancer samples
sample_groups <- metadata %>%
  mutate(
    group = case_when(
      str_detect(str_to_lower(source_name_ch1), "tumor") |
        str_detect(str_to_lower(title), "cancer") ~ "Cancer",

      str_detect(str_to_lower(source_name_ch1), "normal") |
        str_detect(str_to_lower(title), "normal") ~ "Normal",

      TRUE ~ NA_character_
    )
  ) %>%
  select(
    sample_id,
    title,
    source_name_ch1,
    `histopathological exam:ch1`,
    `grade:ch1`,
    `age in years:ch1`,
    race = `race:ch1`,
    group
  )

# -----------------------------
# 3. Check if any sample has no group
# -----------------------------
missing_groups <- sample_groups %>%
  filter(is.na(group))

if (nrow(missing_groups) > 0) {
  cat("WARNING: Some samples could not be labeled.\n")
  print(missing_groups)
} else {
  cat("All samples were labeled successfully.\n")
}

# -----------------------------
# 4. Count samples in each group
# -----------------------------
group_counts <- sample_groups %>%
  count(group)

cat("\nSample counts by group:\n")
print(group_counts)

# -----------------------------
# 5. Save clean sample group table
# -----------------------------
write_csv(sample_groups, "data/processed/GSE15852_sample_groups.csv")
write_csv(group_counts, "results/02_sample_group_counts.csv")

cat("\nFiles created:\n")
cat("data/processed/GSE15852_sample_groups.csv\n")
cat("results/02_sample_group_counts.csv\n")
