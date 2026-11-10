
# ============================================
# ICA 3 - PROBE TO GENE ANNOTATION
# ============================================

library(GEOquery)
library(tidyverse)

# 1. Load differential expression results
results <- read.csv(
  "results/tables/differential_expression_results.csv"
)

# 2. Load GPL6244 annotation
annotation <- read.csv(
  "data/GPL6244_annotation.csv",
  check.names = FALSE
)

# 3. Check annotation columns
cat("Annotation column names:\n")
print(colnames(annotation))

# 4. Convert probe IDs into character format
results$Probe_ID <- as.character(results$Probe_ID)
annotation$ID <- as.character(annotation$ID)

# 5. Find our significant probe
significant_probe <- results %>%
  filter(adj.P.Val < 0.05)

cat("\nSignificant probe:\n")
print(significant_probe)

# 6. Match probe IDs with annotation
annotated_results <- results %>%
  left_join(
    annotation,
    by = c("Probe_ID" = "ID")
  )

# 7. Display annotation for the top probe
top_probe_annotation <- annotated_results %>%
  filter(Probe_ID == "8072314")

cat("\nTop probe annotation:\n")
print(top_probe_annotation)

# 8. Save annotated results
write.csv(
  annotated_results,
  "results/tables/annotated_gene_results.csv",
  row.names = FALSE
)

cat("\nGene annotation matching completed!")
