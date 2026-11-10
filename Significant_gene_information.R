# ============================================
# STEP 1: LOAD ANNOTATED RESULTS
# ============================================

annotated_results <- read.csv(
  "results/tables/annotated_gene_results.csv",
  check.names = FALSE
)

# ============================================
# STEP 2: SELECT SIGNIFICANT PROBES
# ============================================

significant_results <- annotated_results[
  !is.na(annotated_results$adj.P.Val) &
    annotated_results$adj.P.Val < 0.05,
]

# Display significant results
print(significant_results)

# Number of significant probes
cat(
  "Number of significant probes:",
  nrow(significant_results)
)