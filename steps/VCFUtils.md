# VCF Utility Tasks
Utility tasks for VCF file manipulation.  Uses bcftools.

# SortVCFTask
Sorts the input vcf with the option of creating an index file for the final sorted file

## Input Parameters
* File input_vcf - required - the input VCF file
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".sorted" added
* Boolean output_index - optional - If true, an index file will be created for the final sorted file; Defaults to false
* String index_format - optional - The type of format for the index file if output_index is true; Defaults to "csi" with the other option being "tbi"
* String docker_image - required - Docker image to run; must have bcftools in the path
* Int disk_size - optional - disk size to allocate in GB; defaults to 10 plus the input VCF size times 23
* Int mem_size - optional - memory to allocate in GB; defaults to 4 if the input VCF is less than 10GB, 8 otherwise
* Float sort_max_mem - optional - number of GBs for the `--max-mem` parameter; defaults to 0.75 of the memory size
* Int preemptible - optional - the preemptible runtime setting

## Output Parameters
* File output_vcf_gz - the output VCF file
* File output_vcf_idx - the output VCF's index file (if desired)

# FilterVCFWithBEDTask
Filters the VCF to a target region

## Input Parameters
* File input_vcf - required - the input VCF file
* File? input_vcf_idx - optional - the input VCF file index; if specified, filtering is done using `--regions-file`, otherwise `--targets-file` is used
* File input_bed - required - BED file of target regions
* String regions_or_targets - optional - Specifying whether to use targets or regions filtering options; Defaults to an empty string and the WDL will then determine which to filtering option to use
* String target_overlap - optional - the argument for `--regions-overlap` or `--targets-overlap`, defaults to "record"
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".sorted" added
* Boolean output_index - optional - If true, an index file will be created for the final sorted file; Defaults to false
* String index_format - optional - The type of format for the index file if output_index is true; Defaults to "csi" with the other option being "tbi"
* String docker_image - required - Docker image to run; must have bcftools in the path
* Int disk_size - optional - disk size to allocation in GB; defaults to 10 plus the input VCF size times 2 plus the input BED size
* Int preemptible - optional - the preemptible runtime setting
 
## Output Parameters
* File output_vcf_gz - the output VCF file
* File output_vcf_idx - the output VCF's index file (if desired)

# AnnotateVCFTask
Annotates the VCF file

## Input Parameters
* File input_vcf - required - the input VCF file
* File annotations_file - required - VCF, TSV, or BED file of annotations
* File annotations_idx_file - required - index file for annotations
* File headers_file - required - file containing headers to add to the output VCF file
* String column_list - required - column specified for what to add from the annotations file to the VCF file
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".annotated" added
* String docker_image - required - Docker image to run; must have bcftools in the path
* Int disk_size - optional - disk size to allocation in GB; defaults to 10 plus the input VCF size times 2.5 plus the annotations file size
* Int preemptible - optional - the preemptible runtime setting

## Output Parameters
* File output_vcf_gz - the output VCF file

# ConcatVCFsTask
Concatenates a list of VCF files

## Input Parameters
* Array[File] input_vcfs - required - the input VCF files
* Boolean sorted - optional - If true, sort after concatenating; Defaults to false
* Boolean output_index - optional - If true, an index file will be created for the final sorted file; Defaults to false
* String index_format - optional - The type of format for the index file if output_index is true; Defaults to "csi" with the other option being "tbi"
* String output_basename - optional - the basename that will be used to name the output file, defaults to the first input
  VCF name with file extensions removed and ".concat" added
* String docker_image - required - Docker image to run; must have bcftools in the path
* Int disk_size - optional - disk size to allocation in GB; defaults to 10 plus the input VCFs size times 30
* Int mem_size - optional - memory to allocate in GB; defaults to 4 if the input VCFs is less than 10GB, 8 otherwise
* Float sort_max_mem - optional - number of GBs for the `--max-mem` parameter; defaults to 0.75 of the memory size
* Int preemptible - optional - the preemptible runtime setting

## Output Parameters
* File output_vcf_gz - the output VCF file
* File output_vcf_idx - the output VCF's index file (if desired)

# ExtractSamplesFromVCFTask
Extracts one or more samples from a joint VCF file

## Input Parameters
* File input_vcf - required - the input VCF file
* Array[String] sample_ids - required - the list of sample ids to extract
* String min_ac - optional - the minimum allele count option, defaults to "1"
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".subset" added
* String docker_image - required - Docker image to run; must have bcftools in the path
* Int disk_size - optional - disk size to allocation in GB; defaults to 10 plus the input VCF size times 2
* Int preemptible - optional - the preemptible runtime setting

## Output Parameters
* File output_vcf_gz - the output VCF file

# ConvertBCFTask
Converts a BCF file to either vcf.gz or vcf format

## Input Parameters
* File input_bcf - required - the input BCF file
* Boolean output_index - optional - whether or not to create an index for the output VCF; Default is false
* String index_format - optional - format for the output VCF index file (options are "csi" or "tbi"); Default is "csi"
* String output_basename - optional - the basename that will be used to name the output file; defaults to the input BCF's name
  VCF name with file extensions removed and ".subset" added
* String docker_image - required - Docker image to run; must have bcftools in the path
* Int disk_size - optional - disk size to allocate in GB; defaults to 10 plus the input BCF times 2.5
* Int preemptible - optional - the preemptible runtime setting

## Output Parameters
* File output_vcf_gz - the output VCF file
* File output_vcf_idx - the output VCF's index file (if desired)

# MakeCollectiveVCFTask

## Input Parameters
* File input_vcf - required - the input joint VCF file
* String output_basename - optional - the basename that will be used to name the output file; defaults to the input VCF's name plus ".collective"
* String collective_sample_name - optional - name for fake sample; defaults to "CollectiveSample"
* String collective_gt_call - optional - genotype of fake sample, either "0/0", "0/1", "1/1", "./."; defaults to "0/0"
* String docker_image - required - Docker image to run; must be a version of a python docker image
* Int disk_size - optional - disk size to allocate in GB; defaults to 10 plus the input VCF times 1.5
* Int preemptible - optional - the preemptible runtime setting

## Output Parameters
* File output_vcf_gz - the output collective sample VCF file