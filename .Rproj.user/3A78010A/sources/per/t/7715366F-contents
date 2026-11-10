# ============================================
# STEP 1: EXTRACT ENSEMBL GENE ID
# ============================================

ensembl_id <- "ENSG00000268812"

# ============================================
# STEP 2: CHECK AVAILABLE ID TYPES
# ============================================

keytypes(org.Hs.eg.db)

# ============================================
# STEP 3: MAP ENSEMBL TO ENTREZ ID
# ============================================

gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_id,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

# ============================================
# STEP 4: DISPLAY MAPPING RESULT
# ============================================

print(gene_mapping)