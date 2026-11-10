# ============================================
# STEP 1: CHECK GO ANNOTATIONS
# ============================================

go_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "91370",
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GO", "ONTOLOGY", "EVIDENCE")
)

# ============================================
# STEP 2: DISPLAY RESULTS
# ============================================

print(go_annotations)

# ============================================
# STEP 3: COUNT AVAILABLE GO ANNOTATIONS
# ============================================

cat(
  "\nNumber of GO annotations:",
  sum(!is.na(go_annotations$GO))
)