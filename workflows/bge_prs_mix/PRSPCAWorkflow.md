# PRS PCA Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| File | input_vcf | Yes | Joint or single-sample VCF | |
| File | pruning_sites_for_pca | No | Pruning sites used to perform PCA if pca_projections is not provided | |
| String | chr_encoding | No | Chromosome encoding string; used to perform PCA if pca_projections is not provided | |
| String | interaction_docker_image | No | Docker image for use of Python | us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | flash_pca_docker_image | No | Docker image for running flask PCA | us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f |
| String | tidyverse_docker_image | No | Docker image for R tidyverse package | rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1 |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| File | pc_projection | If a PC projections file is not provided | PCA array for adjusting PRS scores |
| File | pc_plot | If a PC projections file is not provided | PCA plot from PCA projection |
