# 04_survival_setup.R — P2 Steps 0-4: build the survival-ready dataset

library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(readxl)
library(survival)
library(caret)
library(TCGAbiolinks)

if (!exists("dds")) dds <- readRDS("dds.rds")
if (!exists("vsd")) vsd <- vst(dds, blind = FALSE)
sig_df <- read.csv("results/deseq2_sig_genes.csv")

# ---- Step 0: survival-ready expression matrix (VST, tumor-only, patient-keyed) ----
tumor_cols <- colData(dds)$sample_type == "Primary Tumor"
expr_tumor <- assay(vsd)[, tumor_cols]
colnames(expr_tumor) <- substr(colnames(expr_tumor), 1, 12)
expr_matrix <- expr_tumor[, !duplicated(colnames(expr_tumor))]
stopifnot(!any(duplicated(colnames(expr_matrix))))
cat(sprintf("Tumor-only expression matrix: %d genes x %d patients\n", nrow(expr_matrix), ncol(expr_matrix)))

# ---- Step 1: survival endpoint from TCGA-CDR ----
cdr <- readxl::read_excel("TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1)
cdr <- cdr[cdr$type == "LUAD", c("bcr_patient_barcode","OS","OS.time","DSS","DSS.time")]
cat(sprintf("CDR rows for LUAD: %d\n", nrow(cdr)))

# ---- Step 2: clinical integration (age -> years; stage cleaned) ----
clinical <- GDCquery_clinic("TCGA-LUAD", type = "clinical")

# GDC's patient-ID column name has drifted across API/package versions —
# detect which one is actually present, fail loudly if neither is.
id_col <- intersect(c("bcr_patient_barcode", "submitter_id"), names(clinical))
stopifnot(length(id_col) >= 1)
clinical$bcr_patient_barcode <- substr(clinical[[id_col[1]]], 1, 12)
stopifnot("ajcc_pathologic_stage" %in% names(clinical))

cdr$bcr_patient_barcode <- substr(cdr$bcr_patient_barcode, 1, 12)

merged <- inner_join(clinical, cdr, by = "bcr_patient_barcode") |>
  filter(!is.na(OS.time), OS.time > 0) |>
  distinct(bcr_patient_barcode, .keep_all = TRUE)

merged$age_years <- if ("age_at_index" %in% names(merged)) {
  as.numeric(merged$age_at_index)
} else {
  as.numeric(merged$age_at_diagnosis) / 365.25
}

merged <- merged[!is.na(merged$ajcc_pathologic_stage) &
                   !merged$ajcc_pathologic_stage %in% c("", "Not Reported", "Stage X"), ]
merged$ajcc_pathologic_stage <- droplevels(factor(merged$ajcc_pathologic_stage))
cat(sprintf("Merged survival cohort: %d patients\n", nrow(merged)))

# ---- Step 3: feature reduction + alignment (unsupervised prefilter + HARD GATE) ----
candidate_genes <- sig_df |> arrange(padj) |> slice_head(n = 500) |> pull(gene)
candidate_genes <- intersect(candidate_genes, rownames(expr_matrix))

common_ids   <- intersect(colnames(expr_matrix), merged$bcr_patient_barcode)
expr_aligned <- expr_matrix[candidate_genes, common_ids, drop = FALSE]
surv_aligned <- merged[match(common_ids, merged$bcr_patient_barcode), ]

stopifnot(identical(colnames(expr_aligned), surv_aligned$bcr_patient_barcode))
cat(sprintf("\n===== ALIGNMENT CHECK =====\nAligned patients: %d | Candidate genes: %d\n",
            ncol(expr_aligned), length(candidate_genes)))

# ---- Step 4: train/test split (stratified, guarded) ----
set.seed(42)
train_idx  <- as.integer(createDataPartition(as.factor(surv_aligned$OS), p = 0.7, list = FALSE))
train_expr <- expr_aligned[,  train_idx, drop = FALSE]
test_expr  <- expr_aligned[, -train_idx, drop = FALSE]
train_surv <- surv_aligned[ train_idx, ]
test_surv  <- surv_aligned[-train_idx, ]

stopifnot(identical(colnames(train_expr), train_surv$bcr_patient_barcode),
          identical(colnames(test_expr),  test_surv$bcr_patient_barcode))

cat(sprintf("\n===== SPLIT CHECK =====\nTrain N: %d | Test N: %d | Train event rate: %.3f\n",
            nrow(train_surv), nrow(test_surv), mean(train_surv$OS)))

saveRDS(list(train_expr=train_expr, test_expr=test_expr, train_surv=train_surv, test_surv=test_surv,
             expr_matrix=expr_matrix, candidate_genes=candidate_genes),
        "p2_setup.rds")
cat("\nSaved p2_setup.rds\n")