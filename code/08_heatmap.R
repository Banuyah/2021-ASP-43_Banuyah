
# ============================================
# ICA 3 - HEATMAP VISUALIZATION
# ============================================

# Install package if needed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)
library(limma)

# 1. Load expression data
eset <- readRDS("data/GSE38642_expression_set.rds")

expression_data <- exprs(eset)
sample_info <- pData(eset)

# 2. Load differential expression results
results <- read.csv(
  "results/tables/differential_expression_results.csv"
)

# 3. Select top 50 probes
top_probes <- results$Probe_ID[
  order(results$adj.P.Val)
]

top_probes <- head(top_probes, 50)

# 4. Extract expression values
heatmap_data <- expression_data[
  as.character(top_probes),
  ,
  drop = FALSE
]

# 5. Create sample annotations
status <- sample_info[["status:ch1"]]

group <- ifelse(
  grepl("T2D", status, ignore.case = TRUE),
  "Type 2 Diabetes",
  "Non-diabetic"
)

annotation_col <- data.frame(
  Group = factor(group)
)

rownames(annotation_col) <- colnames(heatmap_data)

# 6. Create heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Probes"
)

# 7. Save heatmap
png(
  "results/figures/gene_expression_heatmap.png",
  width = 1400,
  height = 1000,
  res = 150
)

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Gene Expression Heatmap"
)

dev.off()

cat("Heatmap generated successfully!")

