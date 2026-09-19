# Snakefile — TCGA-LUAD reproducible pipeline (P1-P4, incl. P2B)
# Run with: snakemake --cores 4 -p

rule all:
    input:
        "figures/pca_plot.png",
        "results/deseq2_sig_genes.csv",
        "figures/volcano_plot.png",
        "p2_setup.rds",
        "p2_models.rds",
        "p2_ph_resolution.rds",
        "p2_setup_v2.rds",
        "p2_validation.rds",
        "p2_refinement.rds",
        "p3b_multiomics.rds",
        "p4_model_selection.rds",
        "p4_grouped_test.rds",
        "models/coef.json",
        "models/model_info.json",
        "models/baseline_surv.csv",
        "gse31210_raw.rds",
        "p2b_external_validation.rds"

rule deseq2_analysis:
    input:
        "data_raw.rds"
    output:
        "dds.rds",
        "figures/pca_plot.png",
        "figures/cook_distance.png"
    log:
        "pipeline_logs/02_deseq2.log"
    shell:
        "Rscript 02_deseq2_analysis.R > {log} 2>&1"

rule deg_enrichment:
    input:
        "dds.rds"
    output:
        "results/deseq2_sig_genes.csv",
        "figures/volcano_plot.png",
        "figures/heatmap_top50.png"
    log:
        "pipeline_logs/03_enrichment.log"
    shell:
        "Rscript 03_deg_enrichment.R > {log} 2>&1"

rule survival_setup:
    input:
        "dds.rds",
        "results/deseq2_sig_genes.csv",
        "TCGA-CDR-SupplementalTableS1.xlsx"
    output:
        "p2_setup.rds"
    log:
        "pipeline_logs/04_survival_setup.log"
    shell:
        "Rscript 04_survival_setup.R > {log} 2>&1"

rule survival_modeling:
    input:
        "p2_setup.rds"
    output:
        "p2_models.rds"
    log:
        "pipeline_logs/05_survival_modeling.log"
    shell:
        "Rscript 05_survival_modeling.R > {log} 2>&1"

rule ph_resolution:
    input:
        "p2_setup.rds",
        "p2_models.rds"
    output:
        "p2_ph_resolution.rds",
        "p2_setup_v2.rds"
    log:
        "pipeline_logs/06_ph_resolution.log"
    shell:
        "Rscript 06_ph_resolution.R > {log} 2>&1"

rule internal_validation:
    input:
        "p2_setup_v2.rds",
        "p2_models.rds"
    output:
        "p2_validation.rds"
    log:
        "pipeline_logs/07_validation.log"
    shell:
        "Rscript 07_internal_validation.R > {log} 2>&1"

rule signature_refinement:
    input:
        "p2_setup_v2.rds",
        "p2_models.rds"
    output:
        "p2_refinement.rds"
    log:
        "pipeline_logs/08_refinement.log"
    shell:
        "Rscript 08_signature_refinement.R > {log} 2>&1"

rule multiomics_tmb:
    input:
        "p2_setup_v2.rds",
        "p2_ph_resolution.rds"
    output:
        "p3b_multiomics.rds"
    log:
        "pipeline_logs/09_multiomics.log"
    retries: 2
    shell:
        "Rscript 09_multiomics_tmb.R > {log} 2>&1"

rule model_selection:
    input:
        "p2_setup_v2.rds",
        "p2_models.rds",
        "p2_ph_resolution.rds"
    output:
        "p4_model_selection.rds"
    log:
        "pipeline_logs/10_model_selection.log"
    shell:
        "Rscript 10_model_selection.R > {log} 2>&1"

rule grouped_stage_test:
    input:
        "p2_setup_v2.rds",
        "p2_ph_resolution.rds"
    output:
        "p4_grouped_test.rds"
    log:
        "pipeline_logs/11_grouped_stage_test.log"
    shell:
        "Rscript 11_grouped_stage_test.R > {log} 2>&1"

rule export_model:
    input:
        "p2_ph_resolution.rds",
        "p2_setup_v2.rds",
        "p4_grouped_test.rds"
    output:
        "models/coef.json",
        "models/model_info.json",
        "models/baseline_surv.csv"
    log:
        "pipeline_logs/12_export_model.log"
    shell:
        "Rscript 12_export_model.R > {log} 2>&1"

rule external_inspect:
    output:
        "gse31210_raw.rds",
        "gse31210_pdata.rds"
    log:
        "pipeline_logs/13a_external_inspect.log"
    retries: 2
    shell:
        "Rscript 13a_external_inspect.R > {log} 2>&1"

rule external_validation:
    input:
        "gse31210_raw.rds",
        "gse31210_pdata.rds",
        "p4_grouped_test.rds",
        "models/coef.json"
    output:
        "p2b_external_validation.rds"
    log:
        "pipeline_logs/13b_external_validation.log"
    shell:
        "Rscript 13b_external_validation.R > {log} 2>&1"
