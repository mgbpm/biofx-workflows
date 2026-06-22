# CRAM to BAM Step
Converts a CRAM file to a coordinate-sorted, indexed BAM using
the supplied reference.  Implemented as a workflow wrapping a
single task so it can be called either as a standalone workflow
or imported as a task from other workflows.

# Task: ConvertCramToBam
Pipes `samtools view` (CRAM→SAM) into `samtools view` (SAM→BAM),
then indexes the result.

# Input Parameters
* File input_cram - required - the CRAM file to convert
* String sample_name - optional - used to name the output files;
  defaults to the input filename with the .cram extension stripped
  (case-insensitive)
* File ref_fasta - required - reference genome FASTA
* File ref_fasta_index - required - reference genome FASTA index
* File ref_dict - required - reference genome dictionary
* String docker - optional - Docker image containing samtools;
  defaults to
  `us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710`
* Int preemptible - optional - preemptible attempt count;
  defaults to `2`

# Output Parameters
* File output_bam - the converted BAM file
* File output_bai - the BAM index file
