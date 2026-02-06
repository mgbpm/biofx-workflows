# RunPRS Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | query_vcf  | Yes | A gz-compressed VCF file of the samples to be scored | |
| File | adjustment_model_manifest | Required to only score the input query VCF | Description of the computed model's parameters in JSON format |
| Array[File] | variant_weights | Required if no model manifest is provided | 3-column TSV files describing variants and their weights | |
| File | score_weights | Required to make a PRSmix model | Score weights for each PGS ID | |
| File | pca_variants | Required if no model manifest is provided | Text file listing the variants for principal component analysis (PCA), one variant per line | |
| File | ref_source | Required if no model manifest is provided | URL to location of reference VCF shards | |
| File | ref_target | Required if no model manifest is provided | URL to use as temporary work space and to save the generated reference VCF | |
| File | condition_code | Yes | Code for condition/disease related to the model and/or scores | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | prs_docker_image | No | Docker image equipped with PRS scripts | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515" |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |
| String | workspace | Yes | name of the Terra workspace where the workflow will be run | |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | raw_scores | Always | Unadjusted PRS for the samples in the query VCF |
| File | mix_score | If PRSmix score is desired | Weighted average of PRS scores |
| File | adjusted_scores | Always | Adjusted PRS single weights score or mix score for the samples in the query VCF |
| File | output_model_manifest | If an input model manifest is not given | Description of the computed model's parameters in JSON format; (this file can be used as one of the inputs for the ScoreQueryVcf workflow) |
| File | kept_pca_variants | If a model manifest is created | File listing the retained and renamed PCA variants |
| Array[File] | renamed_variant_weights | If a model manifest is created | Weights files with renamed variant IDs |
| File | reference_vcf | If a model manifest is created | Generated refernce VCF |
| Array[File] | renamed_query_vcf | If a model manifest is created| Query VCF with renamed variant IDs |
