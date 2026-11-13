
# ==================================================
# ICA 3 - BIOINFORMATICS
# STEP 3: SAMPLE GROUPING AND PCA ANALYSIS
# ==================================================

library(tidyverse)
library(ggplot2)

# 1. Load saved dataset
eset <- readRDS("data/GSE38642_expression_set.rds")

expression_data <- exprs(eset)
sample_info <- pData(eset)

# 2. Match sample information with expression columns
sample_info <- sample_info[
  match(colnames(expression_data), rownames(sample_info)),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(sample_info), colnames(expression_data))
)

# 3. Extract actual diabetes status from NCBI metadata
status <- trimws(
  as.character(sample_info[["status:ch1"]])
)

# 4. Create two sample groups
sample_info$Group <- case_when(
  grepl("non-diabetic", status, ignore.case = TRUE) ~
    "Non-diabetic",
  
  grepl("T2D", status, ignore.case = TRUE) ~
    "Type 2 Diabetes",
  
  TRUE ~ NA_character_
)

# 5. Convert into factor
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Non-diabetic", "Type 2 Diabetes")
)

# Check for unknown labels
if (any(is.na(sample_info$Group))) {
  stop("Some sample groups could not be identified!")
}

# 6. Display group counts
cat("Sample group distribution:\n")
print(table(sample_info$Group))

# Save sample information
write.csv(
  sample_info,
  "data/processed_sample_info.csv",
  row.names = TRUE
)

# ==================================================
# STEP 4: PRINCIPAL COMPONENT ANALYSIS (PCA)
# ==================================================

# Transpose matrix:
# Rows = samples
# Columns = genes/probes

pca_result <- prcomp(
  t(expression_data),
  center = TRUE,
  scale. = FALSE
)

# Calculate percentage variance
variance <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

# Prepare PCA dataframe
pca_df <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = sample_info$Group
)

# Display PCA information
print(head(pca_df))

cat("PC1 variance:", round(variance[1], 2), "%\n")
cat("PC2 variance:", round(variance[2], 2), "%\n")

# ==================================================
# STEP 5: PCA SCATTER PLOT
# ==================================================

p <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Pancreatic Islet Gene Expression",
    subtitle = "Non-diabetic vs Type 2 Diabetes",
    x = paste0("PC1 (", round(variance[1], 2), "%)"),
    y = paste0("PC2 (", round(variance[2], 2), "%)"),
    color = "Sample Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )

print(p)

# Save PCA graph
ggsave(
  "results/figures/PCA_plot.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)

# Save PCA results
write.csv(
  pca_df,
  "results/tables/PCA_results.csv",
  row.names = FALSE
)

cat("\nPCA analysis completed successfully!")

