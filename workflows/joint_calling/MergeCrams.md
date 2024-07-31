# MergeCrams Workflow
A workflow that takes as input paths to two cram files, localizes them, runs cram to bam, merges the bams with samtools merge and converts the output merged bam back to a cram file. This cram file is indexed and both the cram and crai are returned as workflow outputs.

## Input Parameters
## Requirements/expectations :
## - Two analysis-ready CRAM files for a single sample (as identified in RG:SM)
##
## Outputs :
## - One Merged CRAM file and its index

