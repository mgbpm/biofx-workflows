# Make PRS Model Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_code | Yes | Code for condition/disease | |
| Array[File] | variant_weights | Yes | Array of different PGS variant weight files, with 3 columns: variant ID, effect allele, and score | |
| File | pca_variants | Yes | Variants used in PCA projection | |
| File | reference_vcf | Yes | Reference VCF of population for creating/training model | |
| File | query_file | Yes | VCF or TSV to score | |
| File | score_weights | No | Score weights for each PGS ID; required to make a PRSmix model | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | raw_reference_scores | Always | Raw PRS scores from scoring reference VCF with each variant weights file |
| File | mixed_reference_scores | If an input `score_weights` file is supplied | Raw PRS mix scores from reference VCF |
| Array[File] | adjusted_reference_scores | Always | Adjusted reference scores from a single weights PRS model or PRSmix model |
| File | adjustment_model_manifest | Always | JSON of model file paths for a single weights PRS model or PRSmix model |
