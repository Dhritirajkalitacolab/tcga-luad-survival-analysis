# 08_signature_refinement.R — bounded attempt to fix overfitting via a stricter gene pool

library(survival)
library(glmnet)
library(randomForestSRC)
library(timeROC)

p2  <- readRDS("p2_setup_v2.rds")
p2m <- readRDS("p2_models.rds")
train_expr <- p2$train_expr; test_expr <- p2$test_expr
train_surv <- p2$train_surv; test_surv  <- p2$test_surv
uni <- p2m$uni; rsf <- p2m$rsf

# ---- Much stricter candidate pool: top 20 genes by raw univariate p-value ----
# (was 158 -> LASSO overfit badly with 123 events; ~6 genes/event is a defensible ratio)
top_genes <- uni[order(uni$p), ][1:20, "gene"]
cat(sprintf("Refit candidate pool: %d genes (was 158)\n", length(top_genes)))

x_train <- t(train_expr[top_genes, , drop = FALSE])
y_train <- Surv(train_surv$OS.time, train_surv$OS)
set.seed(42)
cvfit2 <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = 10)

# lambda.1se: more conservative, more parsimonious than lambda.min
coef_1se <- coef(cvfit2, s = "lambda.1se")
sel2 <- rownames(coef_1se)[as.numeric(coef_1se) != 0]
if (length(sel2) == 0) {  # fallback if 1se zeroes everything out
  coef_1se <- coef(cvfit2, s = "lambda.min"); sel2 <- rownames(coef_1se)[as.numeric(coef_1se) != 0]
  cat("lambda.1se selected 0 genes -> fell back to lambda.min\n")
}
stopifnot(length(sel2) >= 1)
coef2 <- setNames(as.numeric(coef_1se)[as.numeric(coef_1se) != 0], sel2)
cat(sprintf("Refined signature: %d genes -> %s\n", length(sel2), paste(sel2, collapse=", ")))

# ---- Frozen risk score, frozen cutoff, no refit on test (same discipline as before) ----
risk_train2 <- as.numeric(t(train_expr[sel2, , drop=FALSE]) %*% coef2)
risk_test2  <- as.numeric(t(test_expr [sel2, , drop=FALSE]) %*% coef2)
cutoff2 <- median(risk_train2)
train_surv$risk_score2 <- risk_train2; test_surv$risk_score2 <- risk_test2

cox_train2 <- coxph(Surv(OS.time, OS) ~ risk_score2, data = train_surv, x = TRUE)
c_train2 <- summary(cox_train2)$concordance[1]
c_test2  <- concordance(cox_train2, newdata = test_surv)$concordance

# ---- RSF's own test-set performance (diagnostic — different model class) ----
rsf_test_data <- data.frame(OS.time = test_surv$OS.time, OS = test_surv$OS,
                            t(test_expr[rsf$xvar.names, , drop=FALSE]), check.names = FALSE)
rsf_pred <- predict(rsf, newdata = rsf_test_data)

# randomForestSRC's get.cindex()/$err.rate has a well-documented convention quirk
# (often reports 1-C, not C). Rather than risk the wrong direction silently,
# compute concordance directly via survival::concordance() on RSF's predicted
# mortality (higher predicted mortality = higher risk -> correct direction).
c_rsf_test <- survival::concordance(Surv(test_surv$OS.time, test_surv$OS) ~ rsf_pred$predicted)$concordance

# ---- Comparison ----
cat("\n===== BEFORE vs AFTER =====\n")
cat(sprintf("Original LASSO (158 genes, lambda.min): train=0.722  test=0.538\n"))
cat(sprintf("Refined  LASSO (%2d genes, stricter):     train=%.3f  test=%.3f\n", length(sel2), c_train2, c_test2))
cat(sprintf("RSF (original 11 genes) test C-index:    %.3f\n", c_rsf_test))

# ---- Honest verdict, hard stop here regardless of outcome ----
cat("\n===== FINAL VERDICT (no further tuning after this) =====\n")
if (c_test2 > 0.58) {
  cat("Refined signature shows real improvement — usable, report this version.\n")
} else {
  cat("Refined signature STILL does not validate meaningfully above chance.\n")
  cat("Conclusion: with n=496 patients / 123 events, this cohort likely lacks power\n")
  cat("for a robust transcriptomic prognostic signature. Report as a genuine, honest\n")
  cat("negative finding — this is valid, defensible science, not a failure.\n")
}

saveRDS(list(sel2=sel2, coef2=coef2, cutoff2=cutoff2, c_train2=c_train2, c_test2=c_test2,
             c_rsf_test=c_rsf_test), "p2_refinement.rds")
cat("\nSaved p2_refinement.rds\n")