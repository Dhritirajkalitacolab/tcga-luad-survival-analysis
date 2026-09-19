# 03_deg_enrichment.R — P1 continued: DEG selection, enrichment, volcano + heatmap

library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(enrichplot)

# Reload if this is a fresh session
if (!exists("dds")) dds <- readRDS("dds.rds")
if (!exists("res"))  res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
if (!exists("vsd"))  vsd <- vst(dds, blind = FALSE)

# ---- Step 1: DEG selection (data-frame form; up/down split) ----
res_df <- as.data.frame(res) |> rownames_to_column("gene") |> filter(!is.na(padj))
sig_df <- res_df |> filter(padj < 0.05, abs(log2FoldChange) > 1)
cat(sprintf("Significant DEGs (padj<0.05, |log2FC|>1): %d (expect 3,000-8,000)\n", nrow(sig_df)))

up_df   <- sig_df |> filter(log2FoldChange >  1)
down_df <- sig_df |> filter(log2FoldChange < -1)
cat(sprintf("Up in tumor: %d | Down in tumor: %d\n", nrow(up_df), nrow(down_df)))

dir.create("results", showWarnings = FALSE)
write.csv(sig_df, "results/deseq2_sig_genes.csv", row.names = FALSE)

# ---- Step 2: functional enrichment — up/down SEPARATE, with universe ----
map_ids <- function(ens) bitr(ens, "ENSEMBL", "ENTREZID", org.Hs.eg.db)$ENTREZID

universe_entrez <- map_ids(res_df$gene)
up_entrez       <- map_ids(up_df$gene)
down_entrez     <- map_ids(down_df$gene)

cat("\nRunning GO enrichment (up)...\n")
ego_up   <- enrichGO(up_entrez,   OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = universe_entrez)
cat("Running GO enrichment (down)...\n")
ego_down <- enrichGO(down_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = universe_entrez)

write.csv(as.data.frame(ego_up),   "results/go_enrichment_up.csv",   row.names = FALSE)
write.csv(as.data.frame(ego_down), "results/go_enrichment_down.csv", row.names = FALSE)
cat(sprintf("GO terms — up: %d | down: %d\n", nrow(as.data.frame(ego_up)), nrow(as.data.frame(ego_down))))

# ---- Step 3: visualization ----
EnhancedVolcano(res_df, lab = res_df$gene, x = "log2FoldChange", y = "padj",
                pCutoff = 0.05, FCcutoff = 1, title = "TCGA-LUAD: Tumor vs Normal")
ggplot2::ggsave("figures/volcano_plot.png", width = 9, height = 8)

top50 <- sig_df |> arrange(padj) |> slice_head(n = 50) |> pull(gene)
mat   <- assay(vsd)[top50, ]
mat   <- mat - rowMeans(mat)

png("figures/heatmap_top50.png", width = 900, height = 1000)
pheatmap::pheatmap(mat, annotation_col = as.data.frame(colData(dds)["sample_type"]),
                   show_colnames = FALSE, main = "Top 50 DEGs")
dev.off()

cat("\n===== P1 COMPLETE =====\n")
cat("Check figures/volcano_plot.png and figures/heatmap_top50.png\n")