# RunPRS Workflow

This orchestration workflow invokes two sub-workflows, which together
compute raw and population-adjusted polygenic risk scores (PRS) for
the samples mentioned in a query VCF file.  The first sub-workflow
computes an adjustment model and the second workflows computes PRS,
both raw and adjusted according to the adjustment model.

## Input Parameters

| Type    | Name          | Req'd | Description | Default Value |
| :---    | :---          | :---  | :---        | :--- |
| File    | weights       | Yes   | a 3-column TSV file describing variants and their weights | |
| File    | pca_variants  | Yes   | text file listing the variants for principal component analysis (PCA), one variant per line | |
| File    | reference_vcf | Yes   | gz-compressed VCF file for reference population | |
| String  | model_name    | Yes   | identifying label for the adjustment model | |
| File    | query_vcf     | Yes   | gz-compressed VCF file of the samples to be scored | |
| String  | query_name    | No    | identifying label for the query dataset | basename of query_vcf (minus the ".vcf.gz" extension) |
| Boolean | norename      | No    | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |

### Details
<code>weights</code>
The format for the weights file is three tab-separated columns, and a
row of column headers.  (The workflow does not use the column headers
at all.)  Each non-header row should consist of a variant identifier
in the form &lt;CHROM&gt;:&lt;POS&gt;:&lt;REF&gt;:&lt;ALT&gt;, the
effect allele, and a weight (a floating point number).  As an example,
below are the first few rows of a valid weights file (formatted for
clarity).

| variant_id    | effect_allele | effect_weight |
| :---          | :---          |          ---: |
| 1:752721:A:G  | A             | -3.117514e-05 |
| 1:760912:C:T  | C             | -6.630655e-05 |
| 1:887560:A:C  | A             | -1.508063e-05 |
| 1:888659:T:C  | T             |  4.050415e-05 |
| 1:998395:A:G  | A             |  1.338596e-04 |
| 1:1030633:G:A | A             |  1.240169e-04 |

<code>pca_variants</code>

The format for the pca_variants file is one variant identifier per
line, in the form &lt;CHROM&gt;:&lt;POS&gt;:&lt;REF&gt;:&lt;ALT&gt;.
Only the sites that are specified in this file and that also occur in
the reference and query VCFs will be used to train the adjustment
model.

## Output Parameters

| Type    | Name                      | When   | Description |
| :---    | :---                      | :---   | :---        |
| File    | adjustment_model_manifest | Always | description of the computed model's parameters in JSON format; (this file can be used as one of the inputs for the ScoreQueryVcf workflow) |
| Boolean | converged                 | Always | whether the training of the adjustment model succeeded |
| File    | raw_reference_scores      | Always | unadjusted PRS for the samples in the reference VCF |
| File    | adjusted_reference_scores | Always | adjusted PRS for the samples in the reference VCF |
| File    | raw_query_scores          | Always | unadjusted PRS for the samples in the query VCF |
| File    | adjusted_query_scores     | Always | adjusted PRS for the samples in the query VCF |
| File    | pc_projection             | Always | table of reference PC coordinates for the query samples (TSV) |
| File    | pc_plot                   | Always | plot of the projection of the query samples' PRS onto the space of the first two reference PC (PNG) |
