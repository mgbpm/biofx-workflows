# PRS PCA Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | Imputed VCF or VCF with genome data | |
| File | adjustment_model_manifest | Yes | Adjustment model manifest file from MakeMixModelWorkflow | |
| File | raw_scores | Yes | Raw PRS scores to adjust | |
| File | output_basename | Yes | Basename to use for output file | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| File | renaming_lookup | No | Mapping file for renaming chromosomes | "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | pc_projection | Always | PCA array for adjusting PRS scores |
| File | pc_plot | Always | PCA plot from PCA projection |
| File | adjusted_scores | Always | Adjusted raw PRS scores based on provided adjustment model and PCA |
