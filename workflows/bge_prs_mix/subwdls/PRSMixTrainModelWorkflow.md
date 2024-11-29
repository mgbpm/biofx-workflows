# PRS Mix Train Model Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_name | Yes | Name of condition/disease | |
| Array[File] | var_weights | Yes | Array of different PGS variant weight files, with 3 columns: variant ID, effect allele, and score | |
| File | scoring_sites | Yes | Sites to use in scoring the VCF | |
| File | population_vcf | Yes | VCF to score; from PerformPopulationPCA WDL; variant IDs much match those in var_weights | |
| Int | scoring_mem | No | Memory usage for scoring the input VCF | 8 |
| File | score_weights | Yes | Score weights for each PGS ID in var_weights | |
| File | population_pcs | Yes | Population PCs file from PerformPopulationPCA WDL | |
| File | python_docker_image | No | Python Docker image | "python:3.9.10" |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | ubuntu_docker_image | No | Ubuntu Docker image | "ubuntu:21.10" |
| String | tidyverse_docker_image | No | R tidyverse Docker image | "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | fitted_params | Always | Trained ancestry model parameters |
| Array[File] | raw_population_scores | Always | Raw PRS scores from scoring population VCF with each variant weights file |
| File | mixed_population_scores | Always | Raw PRS mix scores from population VCF |
| File | adjusted_population_scores | Always | Adjusted population scores from the model |
| Boolean | fit_converged | Always | |