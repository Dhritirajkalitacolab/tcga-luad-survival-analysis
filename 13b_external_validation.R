# 13b_external_validation.R — apply the frozen clinical model to GSE31210 (external validation)

library(survival)
library(dplyr)
library(jsonlite)

gset <- readRDS("gse31210_raw.rds")
pdata <- readRDS("gse31210_pdata.rds")
phr <- readRDS("p2_ph_resolution.rds")
sel <- readRDS("p4_grouped_test.rds")
cox_grouped <- sel$cox_grouped
coef_json <- fromJSON("models/coef.json")

# ---- Step 1: filter to tumor tissue only ----
cat("Tissue distribution:\n"); print(table(pdata$`tissue:ch1`))
gse <- pdata[pdata$`tissue:ch1` == "primary lung tumor", ]
cat(sprintf("\nAfter tumor-only filter: %d samples\n", nrow(gse)))

# ---- Step 2: honor the depositors' own prognosis-analysis exclusion flag ----
excl_col <- "exclude for prognosis analysis due to incomplete resection or adjuvant therapy:ch1"
cat("\nExclusion flag distribution:\n"); print(table(gse[[excl_col]], useNA = "always"))
gse <- gse[!is.na(gse[[excl_col]]) & gse[[excl_col]] == "none", ]
cat(sprintf("After honoring depositors' prognosis-analysis exclusion: %d samples\n", nrow(gse)))

# ---- Step 3: extract + type-convert each field ----
gse$age_years <- as.numeric(gse$`age (years):ch1`)
gse$gender    <- as.character(gse$`gender:ch1`)
gse$OS        <- case_when(gse$`death:ch1` == "dead"  ~ 1,
                           gse$`death:ch1` == "alive" ~ 0, TRUE ~ NA_real_)
gse$OS.time   <- as.numeric(gse$`days before death/censor:ch1`)

cat("\nGender values found:\n"); print(table(gse$gender, useNA = "always"))
cat("\nDeath values found:\n");  print(table(gse$`death:ch1`, useNA = "always"))

# ---- Step 4: stage mapping — print FULL distribution first, gate on complete coverage ----
cat("\nFull pathological stage distribution (BEFORE mapping):\n")
print(table(gse$`pathological stage:ch1`, useNA = "always"))

gse$stage_grouped <- case_when(
  gse$`pathological stage:ch1` %in% c("I","IA","IB")     ~ "I",
  gse$`pathological stage:ch1` %in% c("II","IIA","IIB")   ~ "II",
  gse$`pathological stage:ch1` %in% c("III","IIIA","IIIB") ~ "III",
  gse$`pathological stage:ch1` == "IV"                    ~ "IV",
  TRUE ~ NA_character_
)
unmapped <- gse$`pathological stage:ch1`[is.na(gse$stage_grouped) & !is.na(gse$`pathological stage:ch1`)]
if (length(unmapped) > 0) {
  cat("\n!!! UNMAPPED STAGE VALUES FOUND (need to add to case_when):\n")
  print(table(unmapped))
}
gse$stage_grouped <- factor(gse$stage_grouped, levels = c("I","II","III","IV"))

# ---- Step 5: complete cases, unseen-level guard (same pattern as 09/10/11) ----
gse_valid <- gse[complete.cases(gse[, c("age_years","gender","stage_grouped","OS","OS.time")]) &
                   gse$OS.time > 0, ]
cat(sprintf("\n===== FINAL EXTERNAL VALIDATION SET: %d patients =====\n", nrow(gse_valid)))
cat(sprintf("Events: %d (%.1f%%)\n", sum(gse_valid$OS), 100*mean(gse_valid$OS)))

train_stages <- c("I","II","III","IV")  # confirmed levels from 12_export_model.R
unseen <- setdiff(unique(as.character(gse_valid$stage_grouped)), train_stages)
stopifnot(length(unseen) == 0)  # should be empty — I/II/III/IV is a closed set by construction

# ---- Step 6: apply the FROZEN model — no refit — exactly as done for the TCGA test set ----
gse_valid$gender <- factor(gse_valid$gender, levels = c("female","male"))
c_external <- concordance(cox_grouped, newdata = gse_valid)$concordance

# ---- Step 7: redundant cross-check using the raw exported coefficients (matches serve.py logic) ----
lp_manual <- coef_json$age_years * gse_valid$age_years +
  coef_json$gendermale * (gse_valid$gender == "male") +
  coef_json$stage_groupedII  * (gse_valid$stage_grouped == "II")  +
  coef_json$stage_groupedIII * (gse_valid$stage_grouped == "III") +
  coef_json$stage_groupedIV  * (gse_valid$stage_grouped == "IV")
c_manual <- concordance(Surv(gse_valid$OS.time, gse_valid$OS) ~ lp_manual, reverse = TRUE)$concordance

cat("\n===== EXTERNAL VALIDATION RESULT (GSE31210) =====\n")
cat(sprintf("Via cox_grouped object:      C=%.3f\n", c_external))
cat(sprintf("Via raw exported coef.json:  C=%.3f  (cross-check, should match)\n", c_manual))
cat(sprintf("\nFor comparison — TCGA internal held-out test: C=0.608\n"))

saveRDS(list(gse_valid=gse_valid, c_external=c_external, c_manual=c_manual),
        "p2b_external_validation.rds")
cat("\nSaved p2b_external_validation.rds\n")