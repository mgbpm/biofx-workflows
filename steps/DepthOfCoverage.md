# Depth of Coverage Step
The depth of coverage step generates genome wide, per region of interest and per gene coverage metrics.
The step is implemented as a workflow so it can run sub-tasks in parallel and with different container images.

# Sub-tasks
## WGS Task
Runs GATK DepthOfCoverage with no interval file and coverage thresholds of 8, 15 and 30.
## ROI Task
Runs GATK DepthOfCoverage with a single ROI file and coverage thresholds of 8 and 15.
## Gene Task
Runs GATK DepthOfCoverage with an ROI file for each non-overlapping gene cluster and coverage thresholds of 8, 15 and 30.
## Summary Task
Aggregate per gene summaries from the Gene Task and enrich with Entrez IDs.

# Input Parameters
* Boolean run_wgs - optional - whether or not to run the WGS Task, defaults to `true`
* String sample_name - required - name of the sample being processed; used to create output file names
* File ref_fasta - required - reference genome FASTA file
* File ref_fasta_index - required - reference genome FASTA index file
* File ref_dict - required - reference genome dictionary file
* File bam - required - sample BAM file
* File bai - required - sample BAM index file
* String output_basename - optional - the basename to use for output file names, defaults to the bam file with the extension removed
* File? roi_all_bed - optional - BED file of all the regions of interest to analyze; if not specified, the ROI Task is skipped
* Array[RoiAndRefGeneFilePair]? roi_genes - optional - list of ROI and ref gene file pairs; if not specified, the Gene Task is skipped; see biofx-depth-of-coverage/roitools/README.md for tools to generate the files
    * File roi_bed - required - BED file of regions in the gene cluster
    * File ref_gene - required - companion RefGene information file
    * File? ref_gene_idx - optional - index file for the RefGene information file
* File? gene_names - optional - tab-delimited file of gene information, with symbol as the second column and Entrez ID as the second to last column; if not specified, the Summary Task is skipped
* Int gatk_max_heap_gb - optional - the maximum heap size parameter for GATK; defaults to `31`
* Int gatk_disk_size - optional - the size of the disk (gb) for GATK containers; should be large enough to accommodate all the input files; defaults to `100`
* String cov_docker_image - required - the name:tag of the Docker image for the Summary Task
* String gatk_docker_image - optional - the name:tag of the Docker image for GATK3 tasks; defaults to `broadinstitute/gatk3:3.7-0`

# Output Parameters
* File wgs_sample_summary - sample summary output from the WGS Task
* File wgs_sample_statistics - sample statistics output from the WGS Task
* File roi_sample_interval_summary - sample interval summary from the ROI Task
* File roi_sample_interval_statistics - sample interval statistics from the ROI Task
* File roi_sample_statistics - sample statistics from the ROI Task
* File roi_sample_summary - sample summary from the ROI Task
* File roi_sample_cumulative_coverage_counts - sample cumulative coverage counts from the ROI Task
* File roi_sample_cumulative_coverage_proportions - sample cumulative coverage proportions from the ROI Task
* File mt_summary - mitochondrial interval summary extracted from the ROI Task's sample interval summary output
* File gene_summary - sample gene summary aggregated from Gene Task outputs
* File gene_summary_unknown - unknown entries from the aggregated sample gene summary
* File gene_summary_entrez - aggregated sample gene summary enriched with Entrez Gene IDs

# Docker Image Requirements
* WGS, ROI and Genes Tasks - GATK3
* Summary Task
  * python3
  * biofx-depth-of-coverage code repository