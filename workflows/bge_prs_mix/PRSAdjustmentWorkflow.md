# PRS Adjustment Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_name | Yes | Name of condition/disease | |
| File | var_weights | Required if chromosome encoding is not provided | PGS variant weight files, with 3 columns: variant ID, effect allele, and score | |
| File | fitted_model_params | Yes | Adjustment model parameters fitted to input population | |
| File | prs_mix_raw_score | Yes | PRS Mix scores for all subjects | |
| String | weights_chr_encoding | No | Chromosome encoding string; used to perform PCA if pca_projections is not provided | |
| File | pca_projections | Required if inputs to create PCA are not provided | PC projections file | |
| File | input_vcf | Required if PCA projections are not provided | Imputed VCF or VCF with genome data | |
| File | pc_loadings | Required if PCA projections are not provided | | |
| File | pc_meansd | Required if PCA projections are not provided | | |
| File | population_pcs | Required if PCA projections are not provided | PCA for population, the same population used to create the fitted model parameters | |
| File | pruning_sites_for_pca | Required if PCA projections are not provided | Pruning sites used to perform PCA | |
| File | python_docker_image | No | Python Docker image | "python:3.9.10" |
| String | tidyverse_docker_image | No | Docker image for R tidyverse package | rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1 |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | flash_pca_docker_image | No | Docker image for running flask PCA | us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | pc_projection | If a PC projections file is not provided | PCA array for adjusting PRS scores |
| File | pc_plot | If a PC projections file is not provided | PCA plot from PCA projection |
| File | adjusted_scores | Always | Adjusted PRS Mix score based on population model |
