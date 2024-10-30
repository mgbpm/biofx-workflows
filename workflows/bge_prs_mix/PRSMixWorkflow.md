# PRS Mix Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | Imputed VCF or VCF with genome data | |
| File | input_vcf_idx | Yes | Index for input VCF | |
| Array[File] | var_weights_files | Yes | Variant weights files for PRS calculation | |
| File | score_weights_file | Yes | Weights for PRS scores according to PGS ID for calculating PRS Mix scores | |
| File | population_loadings | Yes | Population vectors from PCA | |
| File | population_meansd | Yes | Population statistics from PCA | |
| File | population_pcs | Yes | Principal components from PCA | |
| File | pruning_sites_for_pca | Yes | Pruning sites for PCA | |
| File | scoring_sites | Yes | Sites used in scoring to create ancestry adjustment model | |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- |
| Array[File] | prs_raw_scores | Always | Raw PRS scores for each sample |
| Array[File] | prs_raw_scores_log | Always | Log of raw PRS score calculation |
| Array[File] | prs_raw_sites_scored | Always | Sites from calculating raw PRS scores |
| File | prs_mix_raw_score | Always | Per sample weighted average of raw PRS scores |
| File | adjusted_scores | Always | Adjusted PRS Mix score based on ancestry adjustment model |
| File | pc_projection | Always | PCA array for adjusting PRS scores |
| File | pc_plot | Always | PCA plot from PCA projection |
