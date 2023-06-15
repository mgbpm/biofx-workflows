# VCF Utility Tasks
Utility tasks for VCF file manipulation.  Uses bcftools.

# SortVCFTask {
Sorts the input

## Input Parameters
* File input_vcf - required - the input VCF file
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".sorted" added
* String docker_image - optional - Docker image to run; must have bcftools in the path; defaults to "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
* Int disk_size - optional - disk size to allocation in GB; defaults to 20

## Output Parameters
* File output_vcf_gz - the output VCF file

# FilterVCFWithBEDTask
Filters the VCF to a target region

## Input Parameters
* File input_vcf - required - the input VCF file
* File input_bed - required - the target region
* Boolean use_targets - optional - if true, use the `--targets-file` option, otherwise use the `--regions-file` option; defaults to false
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".sorted" added
* String docker_image - optional - Docker image to run; must have bcftools in the path; defaults to "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
* Int disk_size - optional - disk size to allocation in GB; defaults to 20
 
## Output Parameters
* File output_vcf_gz - the output VCF file

# AnnotateVCFTask
Annotates the VCF file

## Input Parameters
* File input_vcf - required - the input VCF file
* File annotations_file - required - VCF, TSV, or BED file of annotations
* File headers_file - required - file containing headers to add to the output VCF file
* String column_list - required - column specified for what to add from the annotations file to the VCF file
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".sorted" added
* String docker_image - optional - Docker image to run; must have bcftools in the path; defaults to "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
* Int disk_size - optional - disk size to allocation in GB; defaults to 50

## Output Parameters
* File output_vcf_gz - the output VCF file
