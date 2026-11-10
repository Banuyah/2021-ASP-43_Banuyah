# ============================================
# IMPROVED GSEA DOT PLOT
# ============================================

library(enrichplot)
library(ggplot2)
library(stringr)

# Create dot plot
gsea_dotplot <- dotplot(
  gsea_results,
  showCategory = 15,
  title = "GO Biological Process GSEA in Type 2 Diabetes"
)

# Improve pathway labels and formatting
gsea_dotplot <- gsea_dotplot +
  scale_y_discrete(
    labels = function(x) str_wrap(x, width = 35)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  )

# Display plot
print(gsea_dotplot)

# Save high-resolution image
ggsave(
  "results/figures/GO_GSEA_dotplot_improved.png",
  plot = gsea_dotplot,
  width = 14,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Improved GSEA dot plot saved successfully!")

file.exists(
  "results/figures/GO_GSEA_dotplot_improved.png"
)

# ============================================
# GSEA ENRICHMENT PLOT
# ============================================
# Load packages
library(enrichplot)
library(ggplot2)
library(gridExtra)

# Create oxidative phosphorylation plot
gsea_pathway_plot <- gseaplot2(
  gsea_results,
  geneSetID = "GO:0006119",
  title = "Oxidative Phosphorylation",
  pvalue_table = TRUE
)

# Display plot
print(gsea_pathway_plot)

# Save plot
ggsave(
  "results/figures/oxidative_phosphorylation_GSEA.png",
  plot = gsea_pathway_plot,
  width = 11,
  height = 8,
  dpi = 300
)

cat("Pathway enrichment plot saved successfully!")

file.exists(
  "results/figures/oxidative_phosphorylation_GSEA.png"
)


# ============================================
# STEP 1: FINAL GSEA ANALYSIS
# ============================================

library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(123)

gsea_final <- gseGO(
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

# ============================================
# STEP 2: CONVERT RESULTS
# ============================================

gsea_final_table <- as.data.frame(gsea_final)

# ============================================
# STEP 3: COUNT SIGNIFICANT PATHWAYS
# ============================================

cat(
  "Total pathways:",
  nrow(gsea_final_table),
  "\n"
)

cat(
  "Significant pathways:",
  sum(gsea_final_table$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

# ============================================
# STEP 4: DISPLAY TOP 15 PATHWAYS
# ============================================

top_pathways <- gsea_final_table[
  order(gsea_final_table$p.adjust),
  c("Description", "NES", "pvalue", "p.adjust")
]

print(head(top_pathways, 15))

# ============================================
# STEP 5: SAVE FINAL RESULTS
# ============================================

write.csv(
  gsea_final_table,
  "results/tables/GO_GSEA_final_results.csv",
  row.names = FALSE
)

cat("Final GSEA results saved!")