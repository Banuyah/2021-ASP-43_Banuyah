
# ============================================
# ICA 3 - GENE ANNOTATION
# ============================================

library(GEOquery)
library(tidyverse)

# 1. Load differential expression results
results <- read.csv(
  "results/tables/differential_expression_results.csv"
)

# 2. Download GPL6244 platform information
gpl <- getGEO("GPL6244")

# 3. Extract annotation table
annotation <- Table(gpl)

# 4. Inspect annotation columns
cat("Annotation columns:\n")
print(colnames(annotation))

# 5. Display first 10 records
print(head(annotation, 10))

# 6. Save complete annotation
write.csv(
  annotation,
  "data/GPL6244_annotation.csv",
  row.names = FALSE
)

cat("\nPlatform annotation downloaded!")


# ============================================
# FIND TOP STATISTICALLY SIGNIFICANT PROBE
# ============================================

# Sort by adjusted p-value
top_results <- results %>%
  arrange(adj.P.Val)

# Display top 10
print(
  head(
    top_results,
    10
  )
)

# Extract the best-ranked probe
best_probe <- top_results[1, ]

cat("\nBest-ranked probe:\n")
print(best_probe)

cat("\nProbe ID:\n")
print(best_probe$Probe_ID)

cat("\nLog Fold Change:\n")
print(best_probe$logFC)

cat("\nAdjusted P-value:\n")
print(best_probe$adj.P.Val)
