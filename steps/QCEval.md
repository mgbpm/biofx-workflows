# QCEval Task
Applies QC rules depending on project type and adds the QC value to the INFO field in the VCF

## Input Parameters
* File input_vcf - required - the input VCF file
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
  VCF name with file extensions removed and ".qceval" added
* String project_type - required - the type of project, one of "BGE_DRAGEN_TP_BINNING", "WGS", "WGS_DRAGEN", "WES" or "NONE"
* String docker_image - required - the qceval Docker image to use
* Int disk_size - optional - the disk size to allocate in GB, defaults to 12 times the VCF size plus 10

_NOTE:_ the default disk_size assumes that the input VCF is compressed.  For an uncompressed VCF input
the recommended size would be 3.2 times the input VCF plus 10.

## Output Parameters
* File output_vcf_gz - the output VCF file

## Docker Requirements
* bcftools 1.17
* biofx-qceval code repository
* Python 3 with modules: vcfpy