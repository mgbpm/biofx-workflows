# KING

KING (Kinship-based Inference for GWAS) is a relationship inference tool that estimates kinship coefficients infers IBD segments for all pairwise relationships. Unrelated pairs can be precisely separated from close relatives with no false positives, with accuracy up to 3rd- or 4th-degree (depending on array or WGS) for --related and --ibdseg analyses, and up to 2nd-degree for --kinship analysis.

The KING orchestration workflow estimates kinship coefficients from VCF files. There are four run options that can be used to do so: "ibdseg", "kinship", "related", and "duplicate". At least two samples must be used as input to the workflow. Samples for analysis can either be in multiple single-sample VCFs, a single joint VCF, or multiple joint VCFs. This workflow will run with the "related" option by default. (Note: The "ibdseg" option will fail to produce an output if no IBD segments are found. It can also fail if there are less than 10 samples in the analysis, in which case the "kinship" option is better suited.)

## KING Orchestration Input Parameters

An few example JSON input files for running KING are provided in the `example` folder of this repo. The inputs within the JSON files are dummy paths and are not meant to be used as is. Descriptions of each input are outlined below:

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | input_vcfs | Yes | Joint VCFs for identifying related individuals; VCFs must not have overlapping samples and should be gzipped | |
| Array[File] | input_vcfs_idx | No | Index files for input_vcfs; if only one VCF is used, no index file is required | |
| String | output_basename | Yes | Basename for file outputs | |
| File | input_bed | No | BED file for filtering joint VCFs; Use a BED file to increase efficiency of merging dataset VCFs | |
| String | run_type | No | Type of flag to be used for running KING; Either "ibdseg", "kinship", "related", or "duplicate" | "related" |
| Int | degree | No | The maximum degree of relatedness to include in KING output | 3 |
| String | bcftools_docker_image | No | Docker image for bcftools | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17" |
| String | king_docker_iamge | No | Docker image with KING tools | "uwgac/topmed-master@sha256:0bb7f98d6b9182d4e4a6b82c98c04a244d766707875ddfd8a48005a9f5c5481e" |

## KING Orchestration Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | kin_output | When either the --kinship or --related flag is used | .kin file that contains kinship coefficients of individuals |
| File | kin0_output | When either the --kinship or --related flag is used | Second .kin file that contains kinship coefficients of between-family relationship checking |
| File | seg_output | When the --ibdseg flag is used | .seg file that contains kinship coefficients and inferred relationships of samples |
| File | con_output | When the --duplicate flag is used | .con file that contains only duplicate individuals |

## KING WDL Tasks

The WDL tasks used by the KING Orchestration workflow are contained within the KingTasks.wdl document within this repo. This includes tasks that manipulate VCFs and tasks that will run KING. Below are the inputs and outputs for each task:

### FilterVcfTask Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF or VCF gz to filter | |
| File | input_bed | Yes | BED file containing regions to filter to | |
| String | output_basename | No | Basename for output filtered VCF | Defaults to basename of input VCF |
| Int | addldisk | No | Addition disk space to add to the final runtime disk space in GB | 10 |
| Int | preemptible | No | Number of retries for VM | 1 |

### FilterVcfTask Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | output_vcf_gz | If an input BED is supplied | VCF filtered to regions in the input BED file |
| Int | num_snps | Always | Number of SNPs in the input bed file |
| Int | num_samples | Always | Number of samples in the output VCF |

### MergeVcfsTask Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | input_vcfs | Yes | VCFs with non-overlapping samples to merge into one VCF | |
| Array[File] | input_vcfs_idx | No | Index files corresponding to input VCFs; must be in the same order as the input VCF array | |
| File | input_bed | No | BED file with regions to filter the merged VCF to | |
| String | output_basename | Yes | Basename for output files | |
| String | docker_image | Yes | Docker image for bcftools | |
| Int | addldisk | No | Addition disk space to add to the final runtime disk space in GB | 10 |
| Int | mem_size | No | Memory for runtime | Defaults to 4; If the size of input VCFs is greater than 10, defaults to 8 |
| Int | preemptible | No | Number of retries for VM | 2 |

### MergeVcfsTask Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | merged_vcf | Always | Merged VCF of all input VCFs, filtered to regions in the input BED file if given |
| Int | num_snps | Always | Number of SNPs in the input bed file |
| Int | num_samples | Always | Number of samples in the output VCF |

### Vcf2BedTask Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | input_vcf | Yes | VCF to convert to PLINK BED | |
| String | output_basename | No | Basename for output files | Defaults to basename of input VCF |
| String | docker_image | Yes | Docker image that contains PLINK | |
| Int | addldisk | No | Addition disk space to add to the final runtime disk space in GB | 10 |
| Int | plink_mem | No | Memory to use for PLINK in GB; Actual runtime memory will be twice the size of the input PLINK memory | 4 |
| Int | preemptible | No | Number of retries for VM | 1 |

### Vcf2BedTask Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | bed_file | Always | PLINK BED from VCF |
| File | bim_file | Always | BIM file corresponding to output PLINK BED |
| File | fam_file | Always | FAM file corresponding to output PLINK BED |

## RunKingTask

### Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | bed_file | Yes | PLINK BED file from converting input VCF to BED | |
| File | fam_file | Yes | PLINK FAM file corresponding to input BEB | |
| File | bim_file | Yes | PLINK BIM file corresponding to input BED | |
| Int | degree | Yes | Largest degree of relatedness allowed for KING relationships | |
| String | flag | Yes | Flag to run with KING; either "ibdseg", "related", "kinship" or "duplicate" | |
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

## References

Original KING paper:
Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) Robust relationship inference in genome-wide association studies. Bioinformatics 26(22):2867-2873

KING Tutorial: https://www.kingrelatedness.com/manual.shtml
