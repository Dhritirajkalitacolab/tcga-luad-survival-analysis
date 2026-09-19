# 01_download_tcga.R — download TCGA-LUAD RNA-seq (tumor + normal)

library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query, method = "api", files.per.chunk = 20)

# gate: confirm the query resolved to the expected number of files before assembling
n_expected <- nrow(getResults(query))
cat(sprintf("Query resolved to %d files (expect 599)\n", n_expected))

luad_se <- GDCprepare(query, directory = "GDCdata")   # <- renamed from `data`

# sanity check — the numbers that matter
cat("\n===== DOWNLOAD CHECK =====\n")
print(table(colData(luad_se)$sample_type))
cat(sprintf("\nTotal samples: %d\n", ncol(luad_se)))
stopifnot(ncol(luad_se) == n_expected)   # gate: prepared object matches the query

saveRDS(luad_se, "data_raw.rds")
cat("Saved data_raw.rds\n")