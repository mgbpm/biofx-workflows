# KING Orchestration WDL

This workflow will estimate kinship among samples within the input VCFs. The input VCFs must not have overlapping samples in order for the workflow to merge them and perform kinship testing. Input VCFs can be either joint called or single-sample VCFs.

The KING tool itself has several flags to run different relationship inferences, each using a different algorithm or specifications to estimate kinship. The output file types vary for each flag that can be run with KING. Refer to the [KING manual](https://www.kingrelatedness.com/manual.shtml) for further descriptions on each flag.

KING uses PLINK bed, bim, and fam files for relationship inferences. Therefore, this workflow will merge all input VCFs and convert them to these PLINK files before running KING.

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | input_vcfs | Yes | VCFs for identifying related individuals; VCFs within the array must not have overlapping samples and should be gzipped | |
| Array[File] | input_vcfs_idx | Yes | Index files for input_vcfs | |
| File | input_bed | No | BED file for filtering joint VCFs; Use a BED file to increase efficiency of merging dataset VCFs | |
| String | output_basename | Yes | Basename for file outputs | |
| String | run_type | No | Type of flag to be used for running KING; Either "ibdseg", "kinship", "related", or "duplicate" | "related" |
| Int | degree | No | The maximum degree of relatedness to include in KING output | 3 |
| Boolean | missing_to_ref | No | If `true`, all missing variant calls will be converted to reference genotypes (0/0) when merging VCFs | false |
| String | bcftools_docker_image | No | Docker image for bcftools | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | king_docker_iamge | No | Docker image with KING tools | "uwgac/topmed-master@sha256:0bb7f98d6b9182d4e4a6b82c98c04a244d766707875ddfd8a48005a9f5c5481e" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | kin_output | When either the --kinship or --related flag is used | .kin file that contains kinship coefficients of individuals |
| File | kin0_output | When either the --kinship or --related flag is used | Second .kin file that contains kinship coefficients of between-family relationship checking |
| File | seg_output | When the --ibdseg flag is used | .seg file that contains kinship coefficients and inferred relationships of samples |
| File | con_output | When the --duplicate flag is used | .con file that contains only duplicate individuals |
