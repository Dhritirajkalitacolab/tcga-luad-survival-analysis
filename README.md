# TCGA-LUAD Multi-Omic Survival Analysis

A reproducible pipeline that builds and validates a survival model for lung adenocarcinoma (LUAD) patients, using RNA-seq expression, clinical variables, and mutation data from The Cancer Genome Atlas.

## Headline result

**The gene-expression signature did not generalize. The clinical model did — and both results are reported here.**

| Model | Internal C-index | External validation (GSE31210) |
|---|---|---|
| Gene-expression signature (LASSO-Cox, selected from 158 candidate genes) | Train 0.722 → Test **0.538** | Not carried forward — did not generalize |
| Clinical covariate model (age + sex + grouped stage) | Test **0.608** | **0.708** (n = 204, 30 events) |

The LASSO signature looked strong on training data (C=0.722) but collapsed on the held-out test set (C=0.538) — a textbook overfitting signal, given 158 candidate genes against only 123 death events in the whole 496-patient cohort (and fewer in the training split). One bounded refinement was run: shrinking the candidate pool from 158 to the top 20 genes (LASSO kept 5) still failed on the held-out test set (C=0.541), and a Random Survival Forest also failed (test C=0.437). **This is reported as a genuine negative finding**, not hidden or re-tuned until it looked better.

The plain clinical model — age, sex, and a 4-level collapsed AJCC stage — held up. It was deployed, and its coefficients are the ones exported and served by the API in this repo.

## Why the clinical model, not the gene signature, was deployed

Three real engineering decisions drove this, all documented in the scripts:

1. **The original 8-level AJCC staging violated the proportional-hazards assumption** (`cox.zph` p=0.037) — several stage categories had as few as 5 patients. Collapsing to a 4-level grouping (I/II/III/IV) resolved the violation (p=0.141) *and* stabilized the held-out C-index.
2. **A stratified model was compared against the collapsed-covariate model on the same held-out test set** before either was deployed — this wasn't assumed, it was tested (`10_model_selection.R`, `11_grouped_stage_test.R`).
3. **Mutation and tumor-mutation-burden data were tested as an addition** to the validated clinical model via a likelihood-ratio test (`09_multiomics_tmb.R`) — they added no significant prognostic value in this cohort, and that null result is reported rather than omitted.

## Pipeline

```
00_install_packages.R      → environment setup
01_download_tcga.R         → GDC download (599 samples: 540 tumor, 59 normal)
02_deseq2_analysis.R       → DESeq2 differential expression
03_deg_enrichment.R        → DEG selection, GO enrichment, volcano/heatmap
04_survival_setup.R        → survival-ready dataset, train/test split
05_survival_modeling.R     → univariate/multivariate Cox, PH test, LASSO, RSF
06_ph_resolution.R         → diagnoses and resolves the PH violation
07_internal_validation.R   → held-out test: gene signature (0.722 → 0.538)
08_signature_refinement.R  → bounded attempt to fix overfitting (158 → 20 candidate genes)
09_multiomics_tmb.R        → mutation/TMB integration, tested against the clinical model
10_model_selection.R       → held-out comparison: covariate vs. stratified Cox
11_grouped_stage_test.R    → confirms the 4-level grouped model resolves both issues
12_export_model.R          → exports the deployed model's coefficients + baseline curve
13a_external_inspect.R     → downloads and inspects GSE31210 (external cohort)
13b_external_validation.R  → applies the frozen model to GSE31210, cross-checked two ways
```

## External validation details

GSE31210 is an independent Japanese lung adenocarcinoma cohort. Of its 226 primary tumors, 22 are flagged by the depositors for exclusion from prognosis analysis, leaving **204 patients (30 deaths)**. The frozen clinical model was applied without refitting. Its C-index was computed two independent ways — from the R model object and from the exported JSON coefficients that the API uses — and both give **0.708**.

## Reproduce it

```bash
snakemake --cores 4 -p
```

Runs the full DAG end to end: DESeq2 → survival modeling → PH resolution → model export → external validation. Two rules (`multiomics_tmb`, `external_inspect`) depend on live network downloads (GDC, GEO) and are configured with automatic retries, since both failed at least once during real development.

## Serve the model

```bash
docker build -t luad-survival-api .
docker run -p 8000:8000 luad-survival-api
```

`POST /predict` with `age_years`, `gender`, and `stage_grouped` returns a linear predictor, 1/3/5-year survival estimates, and the model's own disclosed limitations. `GET /model-info` surfaces the model's caveats directly from the served API — not just this README.

## Known limitations

- The Stage IV subgroup in the held-out test set is small (n=6, 5 events) — its contribution to the reported C-index should be read with caution.
- Model architecture was selected by comparing three candidates' held-out performance on the same test set; the reported C-index may be marginally optimistic as a result of that selection process.
- **The external cohort contains only stage I and II patients** (IA, IB, II). The external C-index therefore tests the model on early-stage disease only; the stage III and IV coefficients were not tested externally. The external cohort is also small (30 events), so its C-index carries wide uncertainty.
- The gene-expression signature and mutation/TMB profile showed no generalizable prognostic signal beyond the clinical model in this cohort.

## Stack

R (DESeq2, survival, glmnet, randomForestSRC, TCGAbiolinks) · Snakemake · Python (FastAPI) · Docker

**Research/portfolio use only. Not for clinical decision-making.**
