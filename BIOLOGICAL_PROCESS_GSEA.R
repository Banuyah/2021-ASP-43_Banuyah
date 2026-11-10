# ============================================
# ICA 3 - GO BIOLOGICAL PROCESS GSEA
# ============================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# STEP 1: Run GSEA

set.seed(123)

gsea_results <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# STEP 2: Convert results into a data frame

gsea_table <- as.data.frame(gsea_results)

# STEP 3: Display results

cat("Total GO pathways:", nrow(gsea_table), "\n")

print(head(gsea_table, 10))

# STEP 4: Count significant pathways

cat(
  "\nSignificant pathways (adjusted p-value < 0.05):",
  sum(gsea_table$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

# STEP 5: Save results

write.csv(
  gsea_table,
  "results/tables/GO_GSEA_results.csv",
  row.names = FALSE
)

cat("\nGSEA results saved successfully!")