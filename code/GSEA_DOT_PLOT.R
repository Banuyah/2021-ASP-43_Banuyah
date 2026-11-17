# =============================================================
# GO BIOLOGICAL PROCESS GSEA AND ENRICHMENT VISUALISATIONS
# Requires gene_list from ranked_genes_for_GSEA.R
# =============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
})

if (!exists("gene_list", envir = .GlobalEnv) || is.null(gene_list)) {
  stop("gene_list is missing. Run EXTRACT_GENE_IDS.R and ranked_genes_for_GSEA.R first.")
}
if (is.null(names(gene_list)) || anyDuplicated(names(gene_list))) {
  stop("gene_list must have unique Entrez IDs as names.")
}
gene_list <- sort(gene_list[is.finite(gene_list)], decreasing = TRUE)

set.seed(123)
gsea_final <- clusterProfiler::gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE
)
gsea_results <- gsea_final
gsea_final_table <- as.data.frame(gsea_final)

write.csv(
  gsea_final_table,
  "results/tables/GO_GSEA_final_results.csv",
  row.names = FALSE
)
write.csv(
  gsea_final_table[gsea_final_table$p.adjust < 0.05, , drop = FALSE],
  "results/tables/significant_GO_pathways.csv",
  row.names = FALSE
)
write.csv(
  gsea_final_table,
  "results/tables/GO_GSEA_results.csv",
  row.names = FALSE
)

cat("GO BP gene sets tested:", nrow(gsea_final_table), "\n")
cat(
  "GO BP gene sets with adjusted p-value < 0.05:",
  sum(gsea_final_table$p.adjust < 0.05, na.rm = TRUE), "\n"
)

if (nrow(gsea_final_table) > 0) {
  top15 <- head(
    gsea_final_table[order(gsea_final_table$p.adjust),
      c("ID", "Description", "NES", "pvalue", "p.adjust")],
    15
  )
  write.csv(top15, "results/tables/Top_15_GO_pathways.csv", row.names = FALSE)

  p <- enrichplot::dotplot(
    gsea_final,
    showCategory = min(15, nrow(gsea_final_table)),
    x = "NES",
    color = "p.adjust"
  ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    labs(
      title = "Top Enriched GO Biological Processes",
      x = "Normalized Enrichment Score (NES)",
      y = "Biological process",
      color = "Adjusted p-value"
    ) +
    theme_minimal()

  ggsave(
    "results/figures/GO_GSEA_final_dotplot.png",
    plot = p, width = 12, height = 9, dpi = 300
  )

  if ("GO:0006119" %in% gsea_final_table$ID) {
    pathway_plot <- enrichplot::gseaplot2(
      gsea_final, geneSetID = "GO:0006119",
      title = "Oxidative Phosphorylation", pvalue_table = TRUE
    )
    ggsave(
      "results/figures/oxidative_phosphorylation_GSEA.png",
      plot = pathway_plot, width = 11, height = 8, dpi = 300
    )
  }
}

cat("GSEA analysis and plots completed.\n")
