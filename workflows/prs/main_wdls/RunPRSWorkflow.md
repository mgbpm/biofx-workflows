# RunPRS Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | query_vcf  | Yes | A gz-compressed VCF file of the samples to be scored | |
| File | adjustment_model_manifest | No | Description of the computed model's parameters in JSON format; Required to only score the input query VCF |
| Array[File] | variant_weights | No | 3-column TSV files describing variants and their weights; Required if no model manifest is provided | |
| File | score_weights | No | Score weights for each PGS ID; required to make a PRSmix model | |
| File | pca_variants  | No | Text file listing the variants for principal component analysis (PCA), one variant per line; Required if no model manifest is provided | |
| File | reference_vcf | No | A gz-compressed VCF file for reference population; Required if no model manifest is provided | |
| File | condition_code | Yes | Code for condition/disease related to the model and/or scores | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | raw_scores | Always | Unadjusted PRS for the samples in the query VCF |
| File | mix_score | If PRSmix score is desired | Weighted average of PRS scores |
| File | adjusted_scores  | Always | Adjusted PRS single weights score or mix score for the samples in the query VCF |
| File | adjustment_model_manifest | Always | description of the computed model's parameters in JSON format; (this file can be used as one of the inputs for the ScoreQueryVcf workflow) |
