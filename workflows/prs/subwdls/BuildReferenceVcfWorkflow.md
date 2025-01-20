# BuildReferenceVcf Workflow

Build a single reference VCF file that includes only the variants mentioned in the provided weights and pca_variants files.

## Input Parameters

| Type        | Name          | Req'd | Description | Default Value |
| :---        | :---          | :---  | :---        | :---          |
| Array[File] | weights       | Yes   | array of 3-column TSV files describing variants and their weights (the content and format of each of these files is as described for the <code>weights</code> input parameter of the RunPRS workflow) | |
| File        | pca_variants  | Yes   | text file listing the variants for principal component analysis (PCA), one variant per line (the content and format of this file is as described for the <code>pca_variants</code> input parameter of the RunPRS workflow) | |
| String      | workspace     | Yes   | name of the Terra workspace where the workflow will be run | |
| String      | shards_source | Yes   | URL to location of reference VCF shards | |
| String      | work_bucket   | Yes   | URL to use as temporary work space and to save the generated reference VCF | |
| Boolean     | nocleanup     | No    | If true, do NOT rename the chromosomes in the chromosome names in the weights and pca_variants files | false |
| Int         | nbatches      | No    | number of parallel jobs for subsetting the reference shards | 500 |

## Output Parameters

| Type        | Name                  | When   | Description |
| :---        | :---                  | :---   | :---        |
| File        | reference_vcf         | Always | Generated refernce VCF |
| File        | reference_tbi         | Always | Index for generated reference VCF |
| File        | regions               | Always | Regions file used to generate reference VCF |
| Array[File] | renamed_weights       | Always | Weights files with renamed variant IDs |
| File        | renamed_pca_variants  | Always | File listing the retained and renamed PCA variants |
