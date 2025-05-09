# MakeAdjustmentModel Workflow

Generates an adjustment model to adjust raw PRS scores.

## Input Parameters

**NB:** The documentation for the RunPRS workflow's
<code>weights</code>, <code>pca_variants</code>, and
<code>reference_vcf</code> input parameters also applies, without
modification, to this workflow's input parameters of the same name;
see the documentation for the RunPRS workflow for more information.

| Type    | Name          | Req'd | Description | Default Value |
| :---    | :---          | :---  | :---        | :--- |
| File    | weights       | Yes   | a 3-column TSV file describing variants and their weights | |
| File    | pca_variants  | Yes   | text file listing the variants for principal component analysis (PCA), one variant per line | |
| File    | reference_vcf | Yes   | gz-compressed VCF file for reference population | |
| File    | query_file    | Yes   | either a gz-compressed VCF file of the samples to be scored, or the list of variants extracted from one such VCF | |
| String  | name          | Yes   | identifying label for the adjustment model | |
| Boolean | norename      | No    | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |

## Output Parameters

| Type    | Name                      | When   | Description |
| :---    | :---                      | :---   | :---        |
| File    | adjustment_model_manifest | Always | description of the computed model's parameters in JSON format; (this file can be used as one of the inputs for the ScoreQueryVcf workflow) |
| Boolean | converged                 | Always | whether the training of the adjustment model succeeded |
| File    | raw_reference_scores      | Always | unadjusted PRS for the samples in the reference VCF |
| File    | adjusted_reference_scores | Always | adjusted PRS for the samples in the reference VCF |
