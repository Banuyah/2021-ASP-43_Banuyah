# ============================================
# STEP 6: SAVE FINAL GSEA RESULTS AND FIGURES
# ============================================

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)

# 1. Create results folders if necessary
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# 2. Extract top 15 pathways
top15 <- gsea_final_table[
  order(gsea_final_table$p.adjust),
  c("ID", "Description", "NES", "pvalue", "p.adjust")
]

top15 <- head(top15, 15)

# Display top 15
print(top15)

# Save top 15 table
write.csv(
  top15,
  "results/tables/Top_15_GO_pathways.csv",
  row.names = FALSE
)

# 3. Create improved GSEA dot plot
p_gsea <- dotplot(
  gsea_final,
  showCategory = 15,
  x = "NES",
  color = "p.adjust"
) +
  scale_y_discrete(
    labels = function(x) str_wrap(x, width = 40)
  ) +
  labs(
    title = "Top 15 Enriched Biological Processes in T2D",
    x = "Normalized Enrichment Score (NES)",
    y = "Biological Process",
    color = "Adjusted P-value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.y = element_text(size = 9)
  )

# Display plot
print(p_gsea)

# Save final plot
ggsave(
  "results/figures/GO_GSEA_final_dotplot.png",
  plot = p_gsea,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Final GSEA table and figure saved successfully!")

file.exists("results/tables/Top_15_GO_pathways.csv")

file.exists("results/figures/GO_GSEA_final_dotplot.png")