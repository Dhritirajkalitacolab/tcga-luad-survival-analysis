# 13a_external_inspect.R — download GSE31210, inspect its actual structure before assuming anything
 
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
library(GEOquery)
 
options(timeout = 600)
 
cat("Downloading GSE31210 (this may take a few minutes)...\n")
gset <- getGEO("GSE31210", GSEMatrix = TRUE, AnnotGPL = FALSE)
gset <- gset[[1]]
cat(sprintf("Downloaded: %d samples, %d probes\n", ncol(gset), nrow(gset)))
 
# ---- Inspect the actual clinical/phenotype structure — DO NOT assume field names ----
pdata <- pData(gset)
cat("\n===== AVAILABLE PHENOTYPE COLUMNS =====\n")
print(colnames(pdata))
 
cat("\n===== FIRST 3 ROWS (all columns) =====\n")
print(head(pdata, 3))
 
# GEO clinical variables are very often buried in free-text characteristics_ch1.* columns
char_cols <- grep("characteristics_ch1", colnames(pdata), value = TRUE)
cat(sprintf("\n===== 'characteristics_ch1' COLUMNS FOUND: %d =====\n", length(char_cols)))
for (cc in char_cols) {
  cat(sprintf("\n--- %s (first 3 unique values) ---\n", cc))
  print(head(unique(pdata[[cc]]), 3))
}
 
saveRDS(gset, "gse31210_raw.rds")
saveRDS(pdata, "gse31210_pdata.rds")
cat("\nSaved gse31210_raw.rds and gse31210_pdata.rds — DO NOT proceed to harmonization until we've reviewed this output.\n")
