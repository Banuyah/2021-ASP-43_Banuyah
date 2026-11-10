
# =====================================================
# ICA 3 - BIOINFORMATICS PROJECT
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# Using limma
# =====================================================

# 1. LOAD LIBRARIES

library(GEOquery)
library(limma)
library(tidyverse)
library(ggplot2)

# 2. LOAD DATASET

eset <- readRDS(
  "data/GSE38642_expression_set.rds"
)

expression_data <- exprs(eset)
sample_info <- pData(eset)

# 3. MATCH SAMPLE INFORMATION

sample_info <- sample_info[
  match(
    colnames(expression_data),
    rownames(sample_info)
  ),
  ,
  drop = FALSE
]

# 4. IDENTIFY SAMPLE GROUPS

status <- sample_info[["status:ch1"]]

group <- ifelse(
  grepl("T2D", status, ignore.case = TRUE),
  "T2D",
  ifelse(
    grepl("non-diabetic", status, ignore.case = TRUE),
    "Control",
    NA
  )
)

group <- factor(
  group,
  levels = c("Control", "T2D")
)

# Check groups
print(table(group))

if (any(is.na(group))) {
  stop("Unknown sample group found!")
}

# 5. CREATE DESIGN MATRIX

design <- model.matrix(~0 + group)

colnames(design) <- c(
  "Control",
  "T2D"
)

print(design)

# 6. CREATE CONTRAST

contrast_matrix <- makeContrasts(
  T2D - Control,
  levels = design
)

print(contrast_matrix)

# 7. FIT LINEAR MODEL

fit <- lmFit(
  expression_data,
  design
)

# 8. APPLY CONTRAST

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

# 9. EMPIRICAL BAYES ANALYSIS

fit2 <- eBayes(fit2)

# 10. EXTRACT DIFFERENTIALLY EXPRESSED GENES

results <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Add probe IDs
results$Probe_ID <- rownames(results)

# Rearrange columns
results <- results %>%
  select(
    Probe_ID,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B
  )

# 11. CLASSIFY SIGNIFICANT PROBES

results$Regulation <- ifelse(
  results$adj.P.Val < 0.05 &
    results$logFC > 1,
  "Upregulated",
  ifelse(
    results$adj.P.Val < 0.05 &
      results$logFC < -1,
    "Downregulated",
    "Not Significant"
  )
)

# 12. COUNT RESULTS

cat("\nTotal probes analysed:",
    nrow(results), "\n")

cat("\nSignificantly upregulated probes:\n")

print(sum(
  results$Regulation == "Upregulated"
))

cat("\nSignificantly downregulated probes:\n")

print(sum(
  results$Regulation == "Downregulated"
))

cat("\nNot significant:\n")

print(sum(
  results$Regulation == "Not Significant"
))

# 13. SAVE RESULTS

write.csv(
  results,
  "results/tables/differential_expression_results.csv",
  row.names = FALSE
)

# Save significant results
significant_results <- results %>%
  filter(adj.P.Val < 0.05)

write.csv(
  significant_results,
  "results/tables/significant_genes.csv",
  row.names = FALSE
)

# 14. DISPLAY TOP 10 RESULTS

cat("\nTop 10 differentially expressed probes:\n")

print(
  head(results, 10)
)

cat("\nDifferential expression analysis completed!")



# ============================================
# IMPROVED VOLCANO PLOT CLASSIFICATION
# ============================================

results$Regulation <- case_when(
  
  results$adj.P.Val < 0.05 &
    results$logFC > 1 ~ "Strongly Upregulated",
  
  results$adj.P.Val < 0.05 &
    results$logFC < -1 ~ "Strongly Downregulated",
  
  results$adj.P.Val < 0.05 ~ "Significant - Modest Change",
  
  TRUE ~ "Not Significant"
)

# Display counts
print(table(results$Regulation))

# Create volcano plot
ggplot(
  results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, 1e-300)),
    color = Regulation
  )
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Gene Expression in T2D",
    subtitle = "Pancreatic Islet Samples",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Classification"
  ) +
  theme_minimal()

ggsave(
  "results/figures/improved_volcano_plot.png",
  width = 10,
  height = 7,
  dpi = 300
)

