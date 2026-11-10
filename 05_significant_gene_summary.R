
# ============================================
# STEP 1: SIGNIFICANT GENE SUMMARY
# ============================================

library(tidyverse)

# Load differential expression results
results <- read.csv(
  "results/tables/differential_expression_results.csv"
)

# 1. Statistically significant probes
significant_p <- results %>%
  filter(adj.P.Val < 0.05)

# 2. Significant probes with positive fold change
upregulated <- results %>%
  filter(
    adj.P.Val < 0.05,
    logFC > 1
  )

# 3. Significant probes with negative fold change
downregulated <- results %>%
  filter(
    adj.P.Val < 0.05,
    logFC < -1
  )

# Display counts
cat("Total probes:", nrow(results), "\n")

cat(
  "Adjusted P-value < 0.05:",
  nrow(significant_p), "\n"
)

cat(
  "Upregulated probes:",
  nrow(upregulated), "\n"
)

cat(
  "Downregulated probes:",
  nrow(downregulated), "\n"
)

# Display top 10 probes
cat("\nTop 10 probes:\n")

print(
  results %>%
    arrange(adj.P.Val) %>%
    select(
      Probe_ID,
      logFC,
      P.Value,
      adj.P.Val
    ) %>%
    head(10)
)

# Save summary
summary_table <- data.frame(
  Category = c(
    "Total probes",
    "Adjusted P-value < 0.05",
    "Upregulated",
    "Downregulated"
  ),
  Count = c(
    nrow(results),
    nrow(significant_p),
    nrow(upregulated),
    nrow(downregulated)
  )
)

write.csv(
  summary_table,
  "results/tables/significance_summary.csv",
  row.names = FALSE
)
