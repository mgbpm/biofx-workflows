# PRS Train Mix Model Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_code | Yes | Code for condition/disease | |
| Array[File] | var_weights | Yes | Array of different PGS variant weight files, with 3 columns: variant ID, effect allele, and score | |
| File | scoring_sites | Yes | Sites to use in scoring the VCF | |
| File | reference_vcf | Yes | VCF to score; from PerformPopulationPCA WDL; variant IDs much match those in var_weights | |
| File | score_weights | Yes | Score weights for each PGS ID in var_weights | |
| Int | scoring_mem | Yes | Memory in gigabytes to use for scoring reference VCF | |
| File | population_pcs | Yes | Population PCs file from PerformPopulationPCA WDL | |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | ubuntu_docker_image | No | Ubuntu Docker image | "ubuntu:21.10" |
| String | tidyverse_docker_image | No | R tidyverse Docker image | "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | fitted_params | Always | Trained ancestry model parameters |
| Array[File] | sites_used_in_scoring | Always | Sites scored when scoring input VCF |
| Array[File] | adjusted_population_scores | Always | Adjusted PRS mix score from ancestry adjustment model |
| Boolean | fit_converged | Always | |
| Array[File] | raw_population_scores | Always | Raw PRS scores from population VCF |
| File | mixed_population_scores | Always | Raw PRS mix scores from population VCF |
