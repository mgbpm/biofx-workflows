# RunPRS Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | query_vcf  | Yes | A gz-compressed VCF file of the samples to be scored | |
| Array[File] | variant_weights | Yes | 3-column TSV files describing variants and their weights | |
| File | score_weights | No | Score weights for each PGS ID; required to make a PRSmix model | |
| File | pca_variants  | Yes   | Text file listing the variants for principal component analysis (PCA), one variant per line | |
| File | reference_vcf | Yes   | A gz-compressed VCF file for reference population | |
| File | condition_code | Yes | Code for condition/disease | |
| Boolean | norename   | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

### Details

 `weights`
The format for the weights file is three tab-separated columns, and a
row of column headers.  (The workflow does not use the column headers
at all.)  Each non-header row should consist of a variant identifier
in the form &lt;CHROM&gt;:&lt;POS&gt;:&lt;REF&gt;:&lt;ALT&gt;, the
effect allele, and a weight (a floating point number).  As an example,
below are the first few rows of a valid weights file (formatted for
clarity).

| variant_id | effect_allele | effect_weight |
| :---    | :---    |    ---: |
| 1:752721:A:G  | A    | -3.117514e-05 |
| 1:760912:C:T  | C    | -6.630655e-05 |
| 1:887560:A:C  | A    | -1.508063e-05 |
| 1:888659:T:C  | T    |  4.050415e-05 |
| 1:998395:A:G  | A    |  1.338596e-04 |
| 1:1030633:G:A | A    |  1.240169e-04 |

 `pca_variants`

The format for the pca_variants file is one variant identifier per
line, in the form &lt;CHROM&gt;:&lt;POS&gt;:&lt;REF&gt;:&lt;ALT&gt;.
Only the sites that are specified in this file and that also occur in
the reference and query VCFs will be used to train the adjustment
model.

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | raw_scores | Always | Unadjusted PRS for the samples in the query VCF |
| File | mix_score | If PRSmix score is desired | Weighted average of PRS scores |
| File | adjusted_scores  | Always | Adjusted PRS single weights score or mix score for the samples in the query VCF |
| File | adjustment_model_manifest | Always | description of the computed model's parameters in JSON format; (this file can be used as one of the inputs for the ScoreQueryVcf workflow) |
