# =============================================================
# COMPLETE REPRODUCIBLE WORKFLOW
# Run from the repository root in a fresh R session.
# =============================================================

required_packages <- c(
  "Biobase", "GEOquery", "limma", "dplyr", "tidyr", "ggplot2",
  "tibble", "stringr", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "enrichplot", "pheatmap"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the missing packages before running: ",
    paste(missing_packages, collapse = ", "),
    "\nSee README.md for installation instructions."
  )
}

if (!file.exists("data/GSE38642_expression_set.rds")) {
  stop("Run this script from the repository root; the expression-set RDS was not found.")
}

# Keep all generated files in the expected folders.
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# 1. Inspect data and create expression-distribution plot.
source("code/02_data_inspection.R", local = .GlobalEnv)

# 2. Group samples and perform PCA.
source("code/03_sample_grouping_PCA.R", local = .GlobalEnv)

# 3. Differential expression and significance summary.
source("code/04_differential_gene_expression.R", local = .GlobalEnv)
source("code/05_significant_gene_summary.R", local = .GlobalEnv)

# 4. Platform annotation and probe-to-gene mapping.
source("code/06_gene_annotation.R", local = .GlobalEnv)
source("code/07_probe_annotation.R", local = .GlobalEnv)

# 5. Extract valid Entrez identifiers and construct a ranked gene list.
source("code/EXTRACT_GENE_IDS.R", local = .GlobalEnv)
source("code/ranked_genes_for_GSEA.R", local = .GlobalEnv)

if (!exists("gene_list", envir = .GlobalEnv) || length(gene_list) < 10) {
  stop("Too few mapped genes were available for a meaningful GSEA run.")
}

# 6. Run GO Biological Process GSEA and create final tables/plots.
source("code/GSEA_DOT_PLOT.R", local = .GlobalEnv)
source("code/final_outcome.R", local = .GlobalEnv)

# 7. Generate the optional heatmap from the top-ranked probes.
source("code/08_heatmap.R", local = .GlobalEnv)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("\nWorkflow finished. Review regenerated files in results/ before submission.\n")
