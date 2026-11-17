# Heatmap of the 50 most statistically significant probes.
suppressPackageStartupMessages({
  library(Biobase)
  library(pheatmap)
})

eset <- readRDS("data/GSE38642_expression_set.rds")
expression_data <- Biobase::exprs(eset)
sample_info <- Biobase::pData(eset)
results <- read.csv("results/tables/differential_expression_results.csv")

sample_info <- sample_info[
  match(colnames(expression_data), rownames(sample_info)),
  , drop = FALSE
]
status <- as.character(sample_info[["status:ch1"]])
group <- ifelse(
  grepl("T2D", status, ignore.case = TRUE),
  "Type 2 Diabetes", "Non-diabetic"
)

ordered_probes <- results$Probe_ID[order(results$adj.P.Val)]
top_probes <- head(ordered_probes, 50)
top_probes <- intersect(as.character(top_probes), rownames(expression_data))
if (length(top_probes) < 2) stop("Not enough probes available for the heatmap.")

heatmap_data <- expression_data[top_probes, , drop = FALSE]
annotation_col <- data.frame(Group = factor(group))
rownames(annotation_col) <- colnames(heatmap_data)

png("results/figures/gene_expression_heatmap.png", width = 1400, height = 1000, res = 150)
pheatmap::pheatmap(
  heatmap_data, scale = "row", annotation_col = annotation_col,
  show_rownames = FALSE, show_colnames = FALSE,
  clustering_method = "complete",
  main = "Top 50 Ranked Differential-Expression Probes"
)
dev.off()
