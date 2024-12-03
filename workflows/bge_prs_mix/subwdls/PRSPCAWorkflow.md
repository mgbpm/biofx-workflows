# PRS PCA Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_name | Yes | Name of condition/disease | |
| File | input_vcf | Yes | Imputed VCF or VCF with genome data | |
| File | adjustment_model_manifest | Yes | Adjustment model manifest file from MakeMixModelWorkflow | |
| File | prs_raw_scores | Yes | Raw PRS scores for all subjects to adjust | |
| File | python_docker_image | No | Python Docker image | "python:3.9.10" |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | flash_pca_docker_image | No | Docker image for running flask PCA | us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f |
| String | tidyverse_docker_image | No | Docker image for R tidyverse package | rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1 |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | pc_projection | If a PC projections file is not provided | PCA array for adjusting PRS scores |
| File | pc_plot | If a PC projections file is not provided | PCA plot from PCA projection |
| File | adjusted_scores | When raw scores are provided | Adjusted raw PRS scores based on provided adjustment model and PCA |
