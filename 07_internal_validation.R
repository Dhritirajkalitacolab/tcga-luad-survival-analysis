# 07_internal_validation.R — P2 Step 12: internal validation on the held-out test set

library(survival)
library(glmnet)
library(timeROC)
library(dplyr)

p2      <- readRDS("p2_setup_v2.rds")
p2m     <- readRDS("p2_models.rds")
train_expr <- p2$train_expr; test_expr <- p2$test_expr
train_surv <- p2$train_surv; test_surv  <- p2$test_surv
lasso_coef     <- p2m$lasso_coef
selected_genes <- p2m$selected_genes

# ---- Build the risk score: FROZEN coefficients, applied (not refit) to test ----
stopifnot(all(selected_genes %in% rownames(train_expr)))
stopifnot(all(selected_genes %in% rownames(test_expr)))

risk_train <- as.numeric(t(train_expr[selected_genes, , drop = FALSE]) %*% lasso_coef)
risk_test  <- as.numeric(t(test_expr [selected_genes, , drop = FALSE]) %*% lasso_coef)

# cutoff frozen on TRAIN only — never recomputed on test
median_cutoff <- median(risk_train)
train_surv$risk_score <- risk_train; test_surv$risk_score <- risk_test
train_surv$risk_group <- ifelse(risk_train > median_cutoff, "High", "Low")
test_surv$risk_group  <- ifelse(risk_test  > median_cutoff, "High", "Low")

cat("===== RISK GROUP SPLIT (test set) =====\n")
print(table(test_surv$risk_group))

# ---- C-index: train vs test (the core generalization check) ----
# IMPORTANT: c_test uses the TRAIN-fitted coefficient applied to test data —
# never refit on test. Refitting risks a sign flip on a small test set (n=148),
# which would silently deflate/invalidate the reported concordance.
cox_score_train <- coxph(Surv(OS.time, OS) ~ risk_score, data = train_surv, x = TRUE)

c_train <- summary(cox_score_train)$concordance[1]
c_test  <- concordance(cox_score_train, newdata = test_surv)$concordance

cat(sprintf("\n===== C-INDEX =====\nTrain: %.3f\nTest:  %.3f\n", c_train, c_test))
cat(sprintf("Drop (train - test): %.3f %s\n", c_train - c_test,
            ifelse(c_train - c_test > 0.1, "(large drop -> watch for overfitting)", "(reasonable)")))

# ---- Time-dependent AUC (only at horizons within observed follow-up) ----
eval_times <- c(365, 1095, 1825)
eval_times <- eval_times[eval_times < max(test_surv$OS.time, na.rm = TRUE)]
cat(sprintf("\nEvaluating AUC at: %s days\n", paste(eval_times, collapse = ", ")))

troc <- timeROC(T = test_surv$OS.time, delta = test_surv$OS, marker = test_surv$risk_score,
                cause = 1, times = eval_times, iid = TRUE)
cat("\n===== TIME-DEPENDENT AUC (test set) =====\n")
print(troc$AUC)

# ---- Kaplan-Meier by risk group, with logrank test ----
km_fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = test_surv)
logrank <- survdiff(Surv(OS.time, OS) ~ risk_group, data = test_surv)
logrank_p <- 1 - pchisq(logrank$chisq, df = 1)
cat(sprintf("\n===== KM LOGRANK (High vs Low risk, test set) =====\nchisq=%.2f  p=%.4f\n",
            logrank$chisq, logrank_p))

saveRDS(list(risk_train=risk_train, risk_test=risk_test, median_cutoff=median_cutoff,
             c_train=c_train, c_test=c_test, troc=troc, km_fit=km_fit, logrank_p=logrank_p),
        "p2_validation.rds")
cat("\n===== VALIDATION COMPLETE — saved p2_validation.rds =====\n")