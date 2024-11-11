# PRS Raw Score Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | condition_name | Yes | Name of condition/disease related to the variant weight files | |
| Array[File] | var_weights | Yes | Array of different PGS variant weight files, with 3 columns: variant ID, effect allele, and score | |
| File | scoring_sites | Yes | Sites to use in scoring the VCF | |
| File | input_vcf | Yes | Joint or single-sample VCF to score | |
| File | python_docker_image | No | Python Docker image | "python:3.9.10" |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[String] | chromosome_encoding | Always | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' |
| Array[File] | prs_raw_scores | Always | Raw scores for each sample in the input VCF |
| Array[File] | prs_raw_scores_log | Always | Log file for finding raw PRS scores |
| Array[File] | prs_sites_scored | Always | List of sites found for PRS scoring |
