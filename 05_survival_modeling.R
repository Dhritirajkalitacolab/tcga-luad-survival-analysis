# 05_survival_modeling.R — P2 Steps 5-11: univariate/multivariate Cox, PH test, LASSO, RSF

library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(dplyr)

p2 <- readRDS("p2_setup.rds")
train_expr <- p2$train_expr; test_expr <- p2$test_expr
train_surv <- p2$train_surv; test_surv  <- p2$test_surv

# ---- Defensive: GDC's sex/gender column name has drifted across schema versions
# (confirmed: this cohort uses "sex_at_birth", not "gender" or "sex").
sex_col <- intersect(c("gender", "sex", "sex_at_birth"), names(train_surv))
stopifnot(length(sex_col) >= 1)
if (sex_col[1] != "gender") {
  train_surv$gender <- train_surv[[sex_col[1]]]
  test_surv$gender  <- test_surv[[sex_col[1]]]
}

# ---- Step 5: univariate Cox (train only, BH-corrected) ----
cat("Running univariate Cox on", nrow(train_expr), "genes...\n")
uni <- do.call(rbind, lapply(rownames(train_expr), function(g) {
  df <- data.frame(time = train_surv$OS.time, status = train_surv$OS,
                   expr = as.numeric(train_expr[g, ]))
  tryCatch({
    s <- summary(coxph(Surv(time, status) ~ expr, data = df))$coefficients
    data.frame(gene = g, p = s[1, 5], HR = s[1, 2])
  }, error = function(e) NULL)
}))
uni <- uni[stats::complete.cases(uni), ]
uni$padj <- p.adjust(uni$p, method = "BH")
sig_uni  <- uni[uni$padj < 0.05, ]
if (nrow(sig_uni) < 30) sig_uni <- uni[order(uni$p), ][seq_len(min(150, nrow(uni))), ]
cat(sprintf("Univariate genes carried forward: %d\n", nrow(sig_uni)))

# ---- Step 6: multivariate Cox (clinical only) ----
train_surv$ajcc_pathologic_stage <- factor(train_surv$ajcc_pathologic_stage)

clin_complete <- complete.cases(train_surv[, c("age_years", "gender", "ajcc_pathologic_stage")])
cat(sprintf("Clinical Cox uses %d/%d patients (%d dropped for missing data)\n",
            sum(clin_complete), nrow(train_surv), sum(!clin_complete)))
cox_clin <- coxph(Surv(OS.time, OS) ~ age_years + gender + ajcc_pathologic_stage, data = train_surv)
cat("\n===== CLINICAL COX SUMMARY =====\n")
print(summary(cox_clin))

# ---- Step 7: PH assumption test (MANDATORY) ----
ph <- cox.zph(cox_clin)
cat("\n===== PH ASSUMPTION TEST =====\n")
print(ph)

# ---- Step 9: LASSO Cox (train only) ----
x_train <- t(train_expr[sig_uni$gene, , drop = FALSE])
y_train <- Surv(train_surv$OS.time, train_surv$OS)
set.seed(42)
cvfit <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = 10)
coef_min <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(coef_min)[as.numeric(coef_min) != 0]
stopifnot(length(selected_genes) >= 1)
lasso_coef <- setNames(as.numeric(coef_min)[as.numeric(coef_min) != 0], selected_genes)
cat(sprintf("\n===== LASSO RESULT =====\nGenes selected: %d\n", length(selected_genes)))
print(selected_genes)

# ---- Step 10: Random Survival Forest ----
cat("\nRunning RSF (this can take a few minutes)...\n")
rsf_data <- data.frame(OS.time = train_surv$OS.time, OS = train_surv$OS,
                       t(train_expr[selected_genes, , drop = FALSE]), check.names = FALSE)
set.seed(42)
rsf <- rfsrc(Surv(OS.time, OS) ~ ., data = rsf_data, ntree = 500, nodesize = 15,
             importance = TRUE, seed = 42)
vimp <- data.frame(Variable = names(rsf$importance), Importance = rsf$importance) |>
  arrange(desc(Importance))
cat("\n===== RSF TOP VARIABLES =====\n")
print(head(vimp, 10))

saveRDS(list(uni=uni, sig_uni=sig_uni, cox_clin=cox_clin, ph=ph, cvfit=cvfit,
             selected_genes=selected_genes, lasso_coef=lasso_coef, rsf=rsf, vimp=vimp),
        "p2_models.rds")
cat("\n===== P2 CORE COMPLETE — saved p2_models.rds =====\n")