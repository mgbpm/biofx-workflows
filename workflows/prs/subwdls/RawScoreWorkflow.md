# PRS Raw Score Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | Joint or single-sample VCF to score | |
| File | adjustment_model_manifest | Yes | Adjustment model manifest file from MakeMixModelWorkflow | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomesInVcf` on `input_vcf` | false |
| File | renaming_lookup | No | Mapping file for renaming chromosomes | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | prs_raw_scores | Always | Raw scores for each sample in the input VCF |
| Array[File] | prs_raw_scores_log | Always | Log file for finding raw PRS scores |
| Array[File] | sites_scored | Always | List of sites found for PRS scoring |
