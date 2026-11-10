
# ============================================
# ICA 3 - BIOINFORMATICS
# STEP 1: DATASET INSPECTION
# ============================================

library(GEOquery)
library(tidyverse)

# Load downloaded dataset
eset <- readRDS("data/GSE38642_expression_set.rds")

# Extract expression matrix
expression_data <- exprs(eset)

# Extract sample information
sample_info <- pData(eset)

# 1. Display dataset dimensions
cat("Number of genes/probes:",
    nrow(expression_data), "\n")

cat("Number of samples:",
    ncol(expression_data), "\n")

# 2. Display first 6 rows
head(expression_data)

# 3. Display sample information
head(sample_info)

# 4. Display metadata column names
colnames(sample_info)

# 5. Check sample descriptions
print(sample_info$title)

print(sample_info$source_name_ch1)

# 6. Check all available sample characteristics
characteristic_columns <- grep(
  "characteristics_ch1",
  colnames(sample_info),
  value = TRUE
)

for (column in characteristic_columns) {
  cat("\nColumn:", column, "\n")
  print(table(sample_info[[column]], useNA = "ifany"))
}

# 7. Check missing values
cat("\nMissing values in expression data:\n")
print(sum(is.na(expression_data)))

cat("\nMissing values in sample information:\n")
print(sum(is.na(sample_info)))

# 8. Check duplicate sample names
cat("\nDuplicate sample names:\n")
print(sum(duplicated(colnames(expression_data))))

# 9. Check expression summary
summary(as.vector(expression_data))


# ============================================
# STEP 2: QUALITY CONTROL VISUALIZATION
# ============================================

library(ggplot2)

# Convert expression matrix to long format
expression_df <- as.data.frame(expression_data)

expression_df$Gene <- rownames(expression_data)

expression_long <- expression_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  )

# Boxplot of expression values
ggplot(
  expression_long,
  aes(x = Sample, y = Expression)
) +
  geom_boxplot(fill = "skyblue") +
  labs(
    title = "Gene Expression Distribution",
    x = "Samples",
    y = "Expression Values"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 5
    )
  )

# Save figure
ggsave(
  "results/figures/expression_boxplot.png",
  width = 12,
  height = 6
)
