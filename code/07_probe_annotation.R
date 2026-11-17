# Join limma probe-level results to GPL6244 platform annotations.
suppressPackageStartupMessages(library(dplyr))

results <- read.csv(
  "results/tables/differential_expression_results.csv",
  stringsAsFactors = FALSE
)
annotation <- read.csv(
  "data/GPL6244_annotation.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("Probe_ID", "adj.P.Val", "logFC") %in% names(results))) {
  stop("Differential-expression table is missing required columns.")
}
if (!"ID" %in% names(annotation)) {
  stop("GPL6244 annotation table has no ID column.")
}

results$Probe_ID <- as.character(results$Probe_ID)
annotation$ID <- as.character(annotation$ID)
annotated_results <- results %>%
  left_join(annotation, by = c("Probe_ID" = "ID"))

write.csv(
  annotated_results,
  "results/tables/annotated_gene_results.csv",
  row.names = FALSE
)

cat("Annotated probes:", sum(!is.na(annotated_results$gene_assignment)), "\n")
cat("Significant annotated probes:\n")
print(annotated_results %>% filter(adj.P.Val < 0.05))
