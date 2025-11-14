# RunPRS Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | query_vcf  | Yes | A gz-compressed VCF file of the samples to be scored | |
| Array[File] | variant_weights | No | A 3-column TSV files describing variants and their weights; Required if an adjustment model is not provided | |
| File | score_weights | No | Score weights for each PGS ID; required for PRSmix scoring | |
| File | pca_variants  | No | Text file listing the variants for principal component analysis (PCA), one variant per line; Required if an adjustment model is not provided | |
| File | reference_vcf | No | A gz-compressed VCF file for reference population; Required if an adjustment model is not provided | |
| Array[File] | adjustment_model_manifest | No | Array of JSON files that describe computed model's parameters for either single weights models or a PRSmix model | |
| String | condition_code | Yes | Code for condition/disease | |
| Boolean | mix_before_adjustment | No | If `true`, mix raw scores before adjusting them | true |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | raw_scores | Always | Unadjusted PRS for the samples in the query VCF |
| Array[File] | mix_score | If `score_weights` file is supplied as input | Weighted average of PRS scores, either before or after adjusting the raw scores |
| Array[File] | adjusted_scores  | Always | Adjusted PRS single weights score or mix score for the samples in the query VCF |
| Array[File] | output_manifests | Always | Array of JSON files that describe computed model's parameters for either single weights models or a PRSmix model |
