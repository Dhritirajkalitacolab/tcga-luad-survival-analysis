# 09_multiomics_tmb.R — P3B: mutation + TMB integration, tested against the VALIDATED clinical model

library(maftools)
library(TCGAbiolinks)
library(survival)
library(dplyr)

p2  <- readRDS("p2_setup_v2.rds")
phr <- readRDS("p2_ph_resolution.rds")
train_surv <- p2$train_surv; test_surv <- p2$test_surv

# ---- Step 1: download + prepare MAF (mutation) data ----
cat("Querying MAF data for TCGA-LUAD...\n")
maf_query <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation")

maf_results <- getResults(maf_query)
cat("Available variant-calling pipelines:\n")
print(table(maf_results$analysis_workflow_type))
cat(sprintf("\n%d total records — ensemble caller, already deduplicated upstream.\n", nrow(maf_results)))

GDCdownload(maf_query)
maf_data <- GDCprepare(maf_query)
stopifnot(nrow(maf_data) > 0)
laml <- read.maf(maf_data)
cat(sprintf("MAF loaded: %d mutation records, %d unique samples\n",
            nrow(maf_data), length(unique(maf_data$Tumor_Sample_Barcode))))

# ---- Step 2: build gene x patient binary mutation matrix ----
mut_mat <- mutCountMatrix(laml)
colnames(mut_mat) <- substr(colnames(mut_mat), 1, 12)
mut_mat <- mut_mat[, !duplicated(colnames(mut_mat)), drop = FALSE]
mut_bin <- (mut_mat > 0) * 1L

drivers <- intersect(c("EGFR", "KRAS", "STK11", "KEAP1", "TP53"), rownames(mut_bin))
stopifnot(length(drivers) >= 1)
cat(sprintf("Driver genes found in MAF: %s\n", paste(drivers, collapse = ", ")))

# ---- Step 3: TMB (tumor mutation burden) ----
nonsyn_classes <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
                    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
                    "Translation_Start_Site", "Nonstop_Mutation")
nonsyn <- maf_data[maf_data$Variant_Classification %in% nonsyn_classes, ]
tmb <- as.data.frame(table(bcr_patient_barcode = substr(nonsyn$Tumor_Sample_Barcode, 1, 12)))
names(tmb) <- c("bcr_patient_barcode", "n_mut")
tmb$TMB <- tmb$n_mut / 38   # exome size ~38 Mb, standard TCGA convention

# ---- Step 4: assemble molecular covariates, align to BOTH train and test cohorts ----
mol <- data.frame(bcr_patient_barcode = colnames(mut_bin),
                  t(mut_bin[drivers, , drop = FALSE]), check.names = FALSE)
names(mol)[-1] <- paste0(drivers, "_mutated")
mol <- merge(mol, tmb[, c("bcr_patient_barcode", "TMB")], by = "bcr_patient_barcode", all.x = TRUE)

align_molecular <- function(surv_df, mol_df) {
  m <- merge(data.frame(bcr_patient_barcode = surv_df$bcr_patient_barcode),
             mol_df, by = "bcr_patient_barcode", all.x = TRUE)
  m <- m[match(surv_df$bcr_patient_barcode, m$bcr_patient_barcode), ]
  stopifnot(identical(m$bcr_patient_barcode, surv_df$bcr_patient_barcode))
  m[is.na(m)] <- 0   # no mutation call -> assume wild-type / TMB=0
  m
}

train_mol <- align_molecular(train_surv, mol)
test_mol  <- align_molecular(test_surv,  mol)
cat(sprintf("\n===== MOLECULAR COVERAGE =====\nTrain: %d/%d patients matched to MAF data\n",
            sum(train_mol$TMB > 0 | rowSums(train_mol[paste0(drivers,"_mutated")]) > 0), nrow(train_surv)))

# ---- Step 5: does adding TMB + driver mutations improve the VALIDATED clinical model? ----
mol_cols <- c(paste0(drivers, "_mutated"), "TMB")
train_combined <- cbind(train_surv, train_mol[, mol_cols])
test_combined  <- cbind(test_surv,  test_mol[,  mol_cols])

# strata() requires test data to contain ONLY stage levels seen during training —
# createDataPartition() doesn't guarantee this. Verify explicitly rather than
# let concordance() fail opaquely downstream.
train_stages <- levels(droplevels(factor(train_combined$ajcc_pathologic_stage)))
test_stage_levels <- unique(as.character(test_combined$ajcc_pathologic_stage))
unseen_in_train <- setdiff(test_stage_levels, train_stages)

if (length(unseen_in_train) > 0) {
  cat(sprintf("Test set contains stage(s) not seen in training: %s\n",
              paste(unseen_in_train, collapse = ", ")))
  cat(sprintf("Dropping %d test patients with unseen stage levels (cannot score them against a stratified model)\n",
              sum(test_combined$ajcc_pathologic_stage %in% unseen_in_train)))
  test_combined <- test_combined[!test_combined$ajcc_pathologic_stage %in% unseen_in_train, ]
}
test_combined$ajcc_pathologic_stage <- factor(test_combined$ajcc_pathologic_stage, levels = train_stages)
train_combined$ajcc_pathologic_stage <- factor(train_combined$ajcc_pathologic_stage, levels = train_stages)

cat(sprintf("Test set for scoring: %d patients (stage-consistent with training)\n", nrow(test_combined)))

cox_baseline <- coxph(Surv(OS.time, OS) ~ age_years + gender + strata(ajcc_pathologic_stage),
                      data = train_combined, x = TRUE)

mol_formula <- as.formula(paste("Surv(OS.time, OS) ~ age_years + gender + strata(ajcc_pathologic_stage) +",
                                paste(mol_cols, collapse = " + ")))
cox_extended <- coxph(mol_formula, data = train_combined, x = TRUE)

c_baseline_test <- concordance(cox_baseline, newdata = test_combined)$concordance
c_extended_test <- concordance(cox_extended, newdata = test_combined)$concordance

# LRT via anova() is only valid if both nested models were fit on the SAME patients.
stopifnot(cox_baseline$n == cox_extended$n)
cat(sprintf("Both models fit on n=%d patients (LRT comparison is valid)\n", cox_baseline$n))

lrt <- anova(cox_baseline, cox_extended)

cat("\n===== CLINICAL vs CLINICAL+MOLECULAR (test C-index) =====\n")
cat(sprintf("Clinical only (validated, from 06):      test C=%.3f\n", c_baseline_test))
cat(sprintf("Clinical + TMB + driver mutations:        test C=%.3f\n", c_extended_test))
cat("\n===== LIKELIHOOD RATIO TEST (train, does molecular data add fit?) =====\n")
print(lrt)

saveRDS(list(mut_bin = mut_bin, tmb = tmb, drivers = drivers,
             cox_baseline = cox_baseline, cox_extended = cox_extended,
             c_baseline_test = c_baseline_test, c_extended_test = c_extended_test, lrt = lrt),
        "p3b_multiomics.rds")
cat("\n===== P3B COMPLETE — saved p3b_multiomics.rds =====\n")