# Utilities Tasks
General purpose tasks

# FailTask
Exits with a non-zero value to force a workflow failure

## Input Parameters
* String error_message - optional - the error message to print to STDERR before exiting

# MakeBatchesTask
Creates batch lists (files with lists of file bucket locations)

## Input Parameters
* Array[File] - input_files - required - the input files of any type
* Int - batch_size - required - Number of input files to include in each batch
* String output_basename - optional - the basename that will be used to name the output file, defaults to the input
* String docker_image - optional - Docker image to run; must have bcftools in the path; defaults to "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
* Int disk_size - optional - disk size to allocation in GB; defaults to 2 times the input VCF size plus the BED file size plus 10
* Int preemptible - optional - default is 1