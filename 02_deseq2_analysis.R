# 02_deseq2_analysis.R — P1: DESeq2 differential expression (TCGA-LUAD)

library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)

# Reload from disk if this is a fresh session (luad_se is a safe name — not a base R function)
if (!exists("luad_se")) luad_se <- readRDS("data_raw.rds")

# ---- Step 1: select the raw counts assay EXPLICITLY ----
stopifnot("unstranded" %in% assayNames(luad_se))
counts <- assay(luad_se, "unstranded")
storage.mode(counts) <- "integer"
cat("Initial samples:", ncol(counts), "| Initial genes:", nrow(counts), "\n")

# ---- Step 2: data integrity checks ----
stopifnot(!any(is.na(colData(luad_se)$sample_type)))
stopifnot(all(sort(unique(colData(luad_se)$sample_type)) ==
                c("Primary Tumor", "Solid Tissue Normal")))

# ---- Step 3: gene ID cleaning (strip versions, collapse PAR_Y duplicates) ----
rn <- sub("\\..*$", "", rownames(counts))
rownames(counts) <- rn
if (any(duplicated(rn))) counts <- rowsum(counts, group = rn)
storage.mode(counts) <- "integer"
stopifnot(!any(duplicated(rownames(counts))))
cat("Genes after ID cleaning:", nrow(counts), "\n")

# ---- Step 4: pre-filter (on the DESeqDataSet, not the raw matrix) ----
coldata <- as.data.frame(colData(luad_se))
coldata$sample_type <- factor(coldata$sample_type,
                              levels = c("Solid Tissue Normal", "Primary Tumor"))
stopifnot(identical(colnames(counts), rownames(coldata)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ sample_type)

keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
cat("Genes after pre-filter:", nrow(dds), "(expect ~25,000-35,000)\n")

# ---- Step 5: run DESeq2 ----
cat("\nRunning DESeq2 — this is the slow step, let it finish...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
cat("\n===== DESeq2 SUMMARY =====\n")
summary(res)

# ---- Step 6: QC plots — NON-NEGOTIABLE GATE ----
vsd <- vst(dds, blind = FALSE)
dir.create("figures", showWarnings = FALSE)

png("figures/pca_plot.png", width = 800, height = 600)
print(plotPCA(vsd, intgroup = "sample_type"))
dev.off()

png("figures/cook_distance.png", width = 800, height = 600)
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2, main = "Cook's distance")
dev.off()

saveRDS(dds, "dds.rds")
cat("\n===== STOP HERE — CHECK THE PCA PLOT BEFORE CONTINUING =====\n")
cat("Open figures/pca_plot.png. Tumor and Normal points MUST form separate clusters.\n")