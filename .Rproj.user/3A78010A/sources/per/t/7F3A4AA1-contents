# ============================================
# STEP 1: LOAD REQUIRED PACKAGES
# ============================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# ============================================
# STEP 2: VERIFY ENTREZ IDs
# ============================================

valid_entrez <- keys(
  org.Hs.eg.db,
  keytype = "ENTREZID"
)

all_results$ENTREZID <- as.character(
  all_results$ENTREZID
)

all_results$ENTREZID[
  !all_results$ENTREZID %in% valid_entrez
] <- NA

cat(
  "Valid Entrez IDs:",
  sum(!is.na(all_results$ENTREZID)),
  "\n"
)

# ============================================
# STEP 3: REMOVE MISSING VALUES
# ============================================

gsea_data <- all_results[
  !is.na(all_results$ENTREZID) &
    !is.na(all_results$t),
]

# ============================================
# STEP 4: REMOVE DUPLICATE GENES
# Keep the probe with the largest absolute t-value
# ============================================

gsea_data <- gsea_data[
  order(abs(gsea_data$t), decreasing = TRUE),
]

gsea_data <- gsea_data[
  !duplicated(gsea_data$ENTREZID),
]

# ============================================
# STEP 5: CREATE RANKED GENE LIST
# ============================================

gene_list <- gsea_data$t

names(gene_list) <- gsea_data$ENTREZID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

# ============================================
# STEP 6: DISPLAY RESULTS
# ============================================

cat(
  "Total ranked genes:",
  length(gene_list),
  "\n"
)

head(gene_list, 10)