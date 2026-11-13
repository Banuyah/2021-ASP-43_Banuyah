gene_list <- gene_list[
  !is.na(gene_list) &
    !is.na(names(gene_list)) &
    names(gene_list) != ""
]

gene_list <- sort(gene_list, decreasing = TRUE)

cat("Ranked genes:", length(gene_list))

# Check whether GSEA returned results
print(gsea_results)

# Convert results into a data frame
gsea_table <- as.data.frame(gsea_results)

# Count identified pathways
cat(
  "\nTotal GO pathways:",
  nrow(gsea_table),
  "\n"
)

# Display first 10 pathways
head(gsea_table, 10)

if (nrow(gsea_table) > 0) {
  
  cat(
    "Significant pathways (adjusted p-value < 0.05):",
    sum(gsea_table$p.adjust < 0.05, na.rm = TRUE),
    "\n"
  )
  
} else {
  cat("No GO pathways were returned.")
}