# ============================================
# STEP 1: LOAD ANNOTATED RESULTS
# ============================================

all_results <- read.csv(
  "results/tables/annotated_gene_results.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ============================================
# STEP 2: EXTRACT ENTREZ IDs
# ============================================

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

library(stringr)

all_results$ENTREZID <- str_extract(
  all_results$gene_assignment,
  "[0-9]+(?=\\s*(///|$))"
)

# ============================================
# STEP 3: COUNT MAPPED PROBES
# ============================================

cat("Total probes:", nrow(all_results), "\n")

cat(
  "Probes with Entrez IDs:",
  sum(!is.na(all_results$ENTREZID)),
  "\n"
)

cat(
  "Probes without Entrez IDs:",
  sum(is.na(all_results$ENTREZID)),
  "\n"
)

cat(
  "Unique Entrez IDs:",
  length(unique(na.omit(all_results$ENTREZID))),
  "\n"
)

# ============================================
# STEP 4: DISPLAY FIRST 10 RESULTS
# ============================================

head(
  all_results[
    , c("Probe_ID", "ENTREZID", "logFC", "t")
  ],
  10
)