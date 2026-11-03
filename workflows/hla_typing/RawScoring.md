# PRS Raw Scoring Workflow

This workflow will generate PRS scores using PLINK2 and an input weights file.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | Sample VCF | |
| File | weights | Yes | Weights file for scoring PRS; must be have columns in this order 1) varaint ID, 2) effect allele, 3) effect weight | |
| String | plink2_docker_image | No | Docker image for PLINK2 | "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | prs_score | Always | Final .sscore file from PLINK2 that contains raw PRS |
| File | sites_scored | Always | Final .sscore.vars file from PLINK2 that contains varaint IDs from scored variants in output .sscore file |
