# Make Mix Model Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_code | Yes | Code for condition/disease | |
| Array[File] | var_weights | Yes | Array of different PGS variant weight files, with 3 columns: variant ID, effect allele, and score | |
| File | pca_variants | Yes | Variants used in PCA projection | |
| File | reference_vcf | Yes | Reference VCF of population for creating/training model | |
| File | query_file | Yes | VCF or TSV to score | |
| File | score_weights | Yes | Score weights for each PGS ID | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |
| String | tidyverse_docker_image | No | Docker image for R tidyverse package | rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1 |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | adjustment_model_manifest | Always | JSON of model file paths |
| Boolean | fit_converged | Always | |
| Array[File] | raw_reference_scores | Always | Raw PRS scores from scoring reference VCF with each variant weights file |
| File | mixed_reference_scores | Always | Raw PRS mix scores from reference VCF |
| File | adjusted_reference_scores | Always | Adjusted reference scores from the model |
