# PRS Raw Score Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_file | Yes | Tar file of condition-/disease-specific variant weights, score weights, and pca files | |
| File | input_vcf | Yes | Joint or single-sample VCF | |
| String | plink_docker_image | No | Docker image for Plink 2 | us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124 |
| String | interaction_docker_image | No | Docker image for use of Python | us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| String | chromosome_encoding | Always | Chromosome encoding based on if mitochondrial variants are represented by 'chrM' or 'chrMT' |
| File | prs_raw_scores | Always | Raw scores for each sample in the input VCF |
| File | prs_raw_scores_log | Always | Log file for finding raw PRS scores |
| File | prs_sites_scored | Always | List of sites found for PRS scoring |
