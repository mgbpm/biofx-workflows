# Probe Coverage Regions Workflow
Exercises the coverage region inputs (BED file and ROI gene
clusters) by running the DepthOfCoverage ROI, Gene, and Summary
tasks against a sample alignment file.  Accepts either a BAM
(with index) or a CRAM; if the input is a CRAM, it is converted
to BAM via the CramToBamWorkflow step before proceeding.
The WGS task is skipped because it does not use the region
inputs under test.

## Input Validation
* The input file extension must be `.bam` or `.cram`
  (case-insensitive); any other extension fails the workflow.
* BAM input requires a `bai` index file.
* CRAM input must not have a `bai` file; providing one is an
  error (it indicates a misconfigured submission).

## Sub-tasks
### ConvertCramToBam (conditional, CRAM only)
Converts the input CRAM to a coordinate-sorted, indexed BAM.
See [CramToBamWorkflow.md](../../steps/CramToBamWorkflow.md).
### DepthOfCoverageWorkflow
Delegates to the DepthOfCoverage step with `run_wgs = false`.
See [DepthOfCoverage.md](../../steps/DepthOfCoverage.md) for
details on the ROI, Gene, and Summary sub-tasks.

## Input Parameters
* File coveragebed - required - BED file of regions of interest
  for the ROI task
* Array[RoiAndRefGeneFilePair] roigenes - required - list of
  ROI and ref gene file pairs for the Gene task
* File genenames - required - tab-delimited gene name file for
  the Summary task
* File bam_or_cram - required - sample BAM or CRAM file
* File? bai - conditional - BAM index; required when
  bam_or_cram is a BAM, must be omitted when it is a CRAM
* File referencefasta - optional - reference genome FASTA;
  defaults to `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta`
* File referenceindex - optional - reference genome FASTA index;
  defaults to `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai`
* File referencedict - optional - reference genome dictionary;
  defaults to `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict`
* String samtools_image - optional - Docker image for
  ConvertCramToBam; defaults to
  `us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710`
* String gatk_image - optional - Docker image for GATK3 tasks;
  defaults to `broadinstitute/gatk3:3.7-0`
* String cov_image - optional - Docker image for the Summary
  task; defaults to
  `us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/coverage:20230630`

## Output Parameters
* File roi_sample_interval_summary - sample interval summary
  from the ROI task
* File roi_sample_interval_statistics - sample interval
  statistics from the ROI task
* File roi_sample_statistics - sample statistics from the
  ROI task
* File roi_sample_summary - sample summary from the ROI task
* File roi_sample_cumulative_coverage_counts - cumulative
  coverage counts from the ROI task
* File roi_sample_cumulative_coverage_proportions - cumulative
  coverage proportions from the ROI task
* File? mt_summary - mitochondrial interval summary extracted
  from the ROI task output
* File gene_summary - aggregated per-gene summary from the
  Gene task
* File gene_summary_unknown - unknown entries from the
  aggregated gene summary
* File gene_summary_entrez - aggregated gene summary enriched
  with Entrez Gene IDs
