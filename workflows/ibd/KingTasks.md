# KING Tasks

## FilterVcfTask

This task will filter the input VCF file to contain only regions within the input BED file. It will then count the number of SNPs and samples in the resulting VCF. If no input BED is given, the task will simply count the number of SNPs and samples in the input VCF.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF or VCF gz to filter | |
| File | input_vcf_idx | Yes | Index file corresponding to input VCF | |
| File | input_bed | No | BED file containing regions for filtering | |
| String | output_basename | No | Basename for output filtered VCF | Defaults to basename of input VCF |
| Int | addldisk | No | Addition disk space to add to the final runtime disk space in GB | 10 |
| Int | preemptible | No | Number of retries for VM | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | output_vcf_gz | If an input BED is supplied | VCF filtered to regions in the input BED file |
| Int | num_snps | Always | Number of SNPs in the input bed file |
| Int | num_samples | Always | Number of samples in the output VCF |

## MergeVcfsTask

This task will merge all the input VCF files into a single VCF. The input VCF files must not have any overlapping samples. If an BED file is supplied, it will simultaneously filter the VCFs to the regions within the BED file. When merging the VCFs, there is an option to convert missing variant calls for any samples to reference calls. Finally, the task will count the number of SNPs and the number of samples in the resulting merged VCF.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | input_vcfs | Yes | VCFs with non-overlapping samples to merge into one VCF | |
| Array[File] | input_vcfs_idx | No | Index files corresponding to input VCFs; must be in the same order as the input VCF array | |
| File | input_bed | No | BED file with regions for filtering | |
| Boolean | missing_to_ref | No | If `true`, all missing variant calls will be converted to reference genotypes (0/0) | false |
| String | output_basename | Yes | Basename for output files | |
| String | docker_image | Yes | Docker image for bcftools | |
| Int | addldisk | No | Additional disk space to add to the final runtime disk space in GB | 10 |
| Int | mem_size | No | Memory for runtime | Defaults to 4; If the size of input VCFs is greater than 10, defaults to 8 |
| Int | preemptible | No | Number of retries for VM | 2 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | merged_vcf | Always | Merged VCF of all input VCFs, filtered to regions in the input BED file if given |
| Int | num_snps | Always | Number of SNPs in the input bed file |
| Int | num_samples | Always | Number of samples in the output VCF |

## Vcf2BedTask

This task will convert a VCF to PLINK bed, bim, and fam files for use with KING.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF to convert to PLINK BED | |
| String | output_basename | No | Basename for output files | Defaults to basename of input VCF |
| String | docker_image | Yes | Docker image that contains PLINK | |
| Int | addldisk | No | Addition disk space to add to the final runtime disk space in GB | 10 |
| Int | plink_mem | No | Memory to use for PLINK in GB; Actual runtime memory will be twice the size of the input PLINK memory | 4 |
| Int | preemptible | No | Number of retries for VM | 1 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | bed_file | Always | PLINK BED from VCF |
| File | bim_file | Always | BIM file corresponding to output PLINK BED |
| File | fam_file | Always | FAM file corresponding to output PLINK BED |

## RunKingTask

This task will run KING, a kinship estimation tool. This tool has several flags to run different relationship inferences, each using a different algorithm or specifications to estimate kinship. The output file types vary for each flag that can be run with KING. Refer to the [KING manual](https://www.kingrelatedness.com/manual.shtml) for further descriptions on each flag.

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | bed_file | Yes | PLINK BED file from converting input VCF to BED | |
| File | fam_file | Yes | PLINK FAM file corresponding to input BEB | |
| File | bim_file | Yes | PLINK BIM file corresponding to input BED | |
| Int | degree | No | Largest degree of relatedness allowed for KING relationships | 3 |
| String | flag | Yes | Flag to run a specified KING algorithm; either "ibdseg", "related", "kinship" or "duplicate" | |
| String | output_basename | Yes | Basename for output files | |
| String | docker_image | Yes | Docker image for running KING | |
| Int | addldisk | No | Addition disk space to add to the final runtime disk space in GB | 10 |
| Int | cpu | No | CPU for runtime | 2 |
| Int | mem_size | No | Memory for runtime | 4 |
| Int | preemptible | No | Number of retries for VM | 2 |

### Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | kin_output | When either the --kinship or --related flag is used | .kin file that contains kinship coefficients of individuals |
| File | kin0_output | When either the --kinship or --related flag is used | Second .kin file that contains kinship coefficients of between-family relationship checking |
| File | seg_output | When the --ibdseg flag is used | .seg file that contains kinship coefficients and inferred relationships of samples |
| File | con_output | When the --duplicate flag is used | .con file that contains only duplicate individuals |
