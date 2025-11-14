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
| File | mix_model_manifest | If an input `score_weights` file is supplied | JSON of model file paths for a PRSmix model |
| Array[File] | raw_model_manifest | If no input `score_weights` file is supplied | Array of JSON files; each JSON contains model file paths a single weights PRS model |
| Array[File] | raw_reference_scores | Always | Raw PRS scores from scoring reference VCF with each variant weights file |
| File | mixed_reference_scores | If an input `score_weights` file is supplied | Raw PRS mix scores from reference VCF |
| File | adjusted_mix_ref_scores | If an input `score_weights` file is supplied | Adjusted reference scores from a PRSmix model |
| Array[File] | adjusted_ref_scores | If no input `score_weights` file is supplied | Adjusted reference scores from a single weights PRS model |
