# Summarise differential-expression results.
suppressPackageStartupMessages(library(dplyr))

results <- read.csv(
  "results/tables/differential_expression_results.csv",
  stringsAsFactors = FALSE
)

stopifnot(all(c("Probe_ID", "logFC", "adj.P.Val", "Regulation") %in% names(results)))

summary_table <- data.frame(
  Category = c(
    "Total probes", "Adjusted P-value < 0.05",
    "Upregulated", "Downregulated", "Not Significant"
  ),
  Count = c(
    nrow(results),
    sum(results$adj.P.Val < 0.05, na.rm = TRUE),
    sum(results$Regulation == "Upregulated", na.rm = TRUE),
    sum(results$Regulation == "Downregulated", na.rm = TRUE),
    sum(results$Regulation == "Not Significant", na.rm = TRUE)
  )
)

write.csv(summary_table, "results/tables/significance_summary.csv", row.names = FALSE)
print(summary_table)
cat("\nTop ten probes by adjusted p-value:\n")
print(head(results[order(results$adj.P.Val), ], 10))
