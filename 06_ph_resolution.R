# 06_ph_resolution.R — Resolving the PH violation on ajcc_pathologic_stage
# (flagged in 05_survival_modeling.R: chisq=14.915, df=7, p=0.037)

library(survival)
library(survminer)
library(dplyr)

p2 <- readRDS("p2_setup.rds")
train_surv <- p2$train_surv
test_surv  <- p2$test_surv
train_surv$ajcc_pathologic_stage <- factor(train_surv$ajcc_pathologic_stage)

# ---- Persistent fix: p2_setup.rds never had "gender" added (05's fix only
# modified its own in-memory copy, never saved back). Fix once, save once,
# so no future script needs to reapply this. ----
sex_col <- intersect(c("gender", "sex", "sex_at_birth"), names(train_surv))
stopifnot(length(sex_col) >= 1)
if (sex_col[1] != "gender") {
  train_surv$gender <- train_surv[[sex_col[1]]]
  test_surv$gender  <- test_surv[[sex_col[1]]]
}
p2$train_surv <- train_surv
p2$test_surv  <- test_surv
saveRDS(p2, "p2_setup_v2.rds")
cat("Saved p2_setup_v2.rds (gender-corrected copy; original p2_setup.rds left untouched)\n")

# ---- Diagnostic: how sparse is each stage cell? (a likely contributor) ----
cat("\n===== EVENTS PER STAGE =====\n")
print(train_surv |> filter(!is.na(gender), !is.na(age_years)) |>
        group_by(ajcc_pathologic_stage) |>
        summarise(n = n(), events = sum(OS), .groups = "drop"))

# ---- Fix 1: stratified Cox — PH-compliant by construction ----
cox_strat <- coxph(Surv(OS.time, OS) ~ age_years + gender + strata(ajcc_pathologic_stage),
                   data = train_surv)
cat("\n===== STRATIFIED COX (primary, PH-safe) =====\n")
print(summary(cox_strat))

ph_strat <- cox.zph(cox_strat)
cat("\n===== PH TEST — STRATIFIED MODEL =====\n")
print(ph_strat)

# ---- Fix 2: AFT model (blueprint Step 8) — keeps a stage effect estimate
# under a modeling framework that does not assume proportional hazards. ----
aft <- survreg(Surv(OS.time, OS) ~ age_years + gender + ajcc_pathologic_stage,
               data = train_surv, dist = "weibull")
cat("\n===== AFT MODEL (Weibull) — stage effect preserved =====\n")
print(summary(aft))

# ---- Diagnostic: does collapsing stage granularity resolve it? ----
train_surv$stage_grouped <- case_when(
  train_surv$ajcc_pathologic_stage %in% c("Stage I","Stage IA","Stage IB")     ~ "I",
  train_surv$ajcc_pathologic_stage %in% c("Stage II","Stage IIA","Stage IIB")  ~ "II",
  train_surv$ajcc_pathologic_stage %in% c("Stage III","Stage IIIA","Stage IIIB") ~ "III",
  train_surv$ajcc_pathologic_stage == "Stage IV"                              ~ "IV",
  TRUE ~ NA_character_
)
stopifnot(sum(is.na(train_surv$stage_grouped)) == 0)
train_surv$stage_grouped <- factor(train_surv$stage_grouped, levels = c("I","II","III","IV"))

cox_grouped <- coxph(Surv(OS.time, OS) ~ age_years + gender + stage_grouped, data = train_surv)
ph_grouped  <- cox.zph(cox_grouped)
cat("\n===== PH TEST — COLLAPSED STAGE (I/II/III/IV) =====\n")
print(ph_grouped)

# ---- Verdict (robust row lookup — doesn't assume exact row-name format) ----
stage_rows <- grep("stage_grouped", rownames(ph_grouped$table), value = TRUE)
stopifnot(length(stage_rows) >= 1)
grouped_p <- max(ph_grouped$table[stage_rows, "p"])

cat("\n===== VERDICT =====\n")
cat(sprintf("Granular stage PH p-value:  0.037 (violated)\n"))
cat(sprintf("Collapsed stage PH p-value: %.3f (%s)\n",
            grouped_p, ifelse(grouped_p > 0.05, "OK — sparsity was a likely driver", "still violated")))
cat("Reporting strategy: use the STRATIFIED Cox model as the primary clinical model\n")
cat("(age/gender HRs valid, PH-compliant by design). Report stage's direction and\n")
cat("magnitude from the AFT model, with the original PH caveat disclosed in the write-up.\n")

saveRDS(list(cox_strat=cox_strat, ph_strat=ph_strat, aft=aft,
             cox_grouped=cox_grouped, ph_grouped=ph_grouped),
        "p2_ph_resolution.rds")
cat("\nSaved p2_ph_resolution.rds\n")