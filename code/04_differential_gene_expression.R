# =============================================================
# ICA 3 - DIFFERENTIAL GENE EXPRESSION ANALYSIS
# Dataset: GSE38642 | Comparison: T2D - Control
# =============================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(dplyr)
  library(ggplot2)
})

# Run this script from the repository root directory.
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

eset <- readRDS("data/GSE38642_expression_set.rds")
expression_data <- Biobase::exprs(eset)
sample_info <- Biobase::pData(eset)

# Align metadata with expression-matrix columns.
sample_info <- sample_info[
  match(colnames(expression_data), rownames(sample_info)),
  , drop = FALSE
]
if (anyNA(rownames(sample_info)) ||
    !identical(rownames(sample_info), colnames(expression_data))) {
  stop("Sample metadata could not be aligned with expression columns.")
}

if (!"status:ch1" %in% colnames(sample_info)) {
  stop("Expected metadata field 'status:ch1' was not found.")
}

status <- trimws(as.character(sample_info[["status:ch1"]]))
group <- dplyr::case_when(
  grepl("non-diabetic", status, ignore.case = TRUE) ~ "Control",
  grepl("T2D", status, ignore.case = TRUE) ~ "T2D",
  TRUE ~ NA_character_
)
if (anyNA(group)) stop("One or more samples have an unrecognised disease status.")
group <- factor(group, levels = c("Control", "T2D"))
print(table(group))

# Fit limma model and estimate T2D minus Control.
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast_matrix <- limma::makeContrasts(T2D - Control, levels = design)
fit <- limma::lmFit(expression_data, design)
fit <- limma::contrasts.fit(fit, contrast_matrix)
fit <- limma::eBayes(fit)

results <- limma::topTable(
  fit, coef = 1, number = Inf,
  adjust.method = "BH", sort.by = "P"
) %>%
  tibble::rownames_to_column("Probe_ID") %>%
  dplyr::select(Probe_ID, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# Direction is assigned only to probes that pass the adjusted p-value cutoff.
# Positive logFC means higher expression in T2D for the T2D - Control contrast.
results <- results %>%
  mutate(
    Regulation = case_when(
      !is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC > 0 ~ "Upregulated",
      !is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

significant_results <- results %>% filter(adj.P.Val < 0.05)
summary_table <- tibble::tibble(
  Category = c(
    "Total probes", "Adjusted P-value < 0.05",
    "Upregulated", "Downregulated", "Not Significant"
  ),
  Count = c(
    nrow(results), nrow(significant_results),
    sum(results$Regulation == "Upregulated"),
    sum(results$Regulation == "Downregulated"),
    sum(results$Regulation == "Not Significant")
  )
)

write.csv(results, "results/tables/differential_expression_results.csv", row.names = FALSE)
write.csv(significant_results, "results/tables/significant_genes.csv", row.names = FALSE)
write.csv(summary_table, "results/tables/significance_summary.csv", row.names = FALSE)

# Volcano plot: show statistical significance and direction without imposing
# an arbitrary |log2 fold-change| > 1 threshold.
volcano_plot <- ggplot(
  results,
  aes(x = logFC, y = -log10(pmax(adj.P.Val, .Machine$double.xmin)), color = Regulation)
) +
  geom_point(alpha = 0.7, size = 1.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Differential Expression in Pancreatic Islets",
    subtitle = "Type 2 Diabetes versus non-diabetic controls",
    x = "Log2 fold change (T2D - Control)",
    y = "-Log10 adjusted p-value",
    color = "Classification"
  ) +
  theme_minimal()

ggplot2::ggsave(
  "results/figures/improved_volcano_plot.png",
  plot = volcano_plot, width = 10, height = 7, dpi = 300
)

cat("\nDifferential-expression analysis completed.\n")
cat("Total probes:", nrow(results), "\n")
cat("Significant probes:", nrow(significant_results), "\n")
print(summary_table)
