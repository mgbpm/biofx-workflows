# PreparePrsMixInputs Workflow

This WDL harmonizes a set of inputs to meet the requirements of the
PRS Mix workflow.  It also generates the smallest reference VCF that
is compatible with the remaining input parameters.

## Input Parameters

| Type        | Name          | Req'd | Description | Default Value |
| :---        | :---          | :---  | :---        | :---          |
| Array[File] | weights_files | Yes   | array of 3-column TSV files describing variants and their weights (the content and format of each of these files is as described for the <code>weights</code> input parameter of the RunPRS workflow) | |
| File        | pca_variants  | Yes   | text file listing the variants for principal component analysis (PCA), one variant per line (the content and format of this file is as described for the <code>pca_variants</code> input parameter of the RunPRS workflow) | |
| String      | workspace     | Yes   | name of the Terra workspace where the workflow will be run | |
| String      | source        | Yes   | URL to location of reference VCF shards | |
| String      | target        | Yes   | URL to use as temporary work space and to save the generated reference VCF | |
| Int         | nbatches      | No    | number of parallel jobs for subsetting the reference shards | 500 |
| Boolean     | resuming      | No    | whether this run is the resumption of an earlier run | false |
| Array[File] | query_vcfs    | Yes   | array of gz-compressed VCF files of the samples to be scored | |

## Output Parameters

| Type        | Name                  | When   | Description |
| :---        | :---                  | :---   | :---        |
| File        | regions               | Always | Regions file used to generate reference VCF |
| File        | kept_pca_variants     | Always | File listing the retained and renamed PCA variants |
| Array[File] | renamed_weights_files | Always | Weights files with renamed variant IDs |
| Array[File] | renamed_query_vcfs    | Always | Query VCFs with renamed variant IDs |
| File        | reference_vcf         | Always | Generated refernce VCF |
| File        | reference_tbi         | Always | Index for generated reference VCF |
