# 00_install_packages.R — install all 17 packages the TCGA-LUAD pipeline needs

cran_pkgs <- c("survival", "survminer", "glmnet", "timeROC", "pec", "cmprsk",
               "rms", "dplyr", "tibble", "jsonlite")

bioc_pkgs <- c("TCGAbiolinks", "DESeq2", "SummarizedExperiment",
               "clusterProfiler", "org.Hs.eg.db", "enrichplot", "EnhancedVolcano")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(setdiff(cran_pkgs, rownames(installed.packages())))
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())), update = FALSE, ask = FALSE)

# verify: every package must load cleanly
all_pkgs <- c(cran_pkgs, bioc_pkgs)
results <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
cat("\n===== INSTALL CHECK =====\n")
print(results)
cat(sprintf("\n%d/%d packages OK\n", sum(results), length(results)))
if (any(!results)) cat("FAILED:", paste(names(results)[!results], collapse=", "), "\n")