# BGWGS Workflow
The BGWGS (bigwigs) Workflow starts with a single CRAM or BAM file and provides variant calling,
filtration, annotation, coverage, pharmacogenetic, risk and reporting capabilities.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | gcp_project_id | Yes | The GCP project to fetch secrets from | |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | Yes | The name of the orchestration utils Docker image for FAST and file movement tasks | |
| String | bcftools_docker_image | Yes | The name of the bcftools Docker image for VCF annotation | |
| String | subject_id | Yes | The subject id associated with the data | |
| String | sample_id | Yes | The sample id associated with the data | |
| String | sample_data_location | Yes | The cloud storage URL where the sample source data is located | |
| Boolean | fetch_cram | No | Whether or not to fetch the CRAM (primarily for testing) | true |
| Boolean | fetch_bam | No | Whether or not to fetch the BAM (primarily for testing) | true |
| Boolean | do_variant_calling | No | Whether or not to generate a VCF by calling variants from BAM or CRAM; if false, the VCF is fetched from the sample_data_location | true |
| String | reference_build | No | The genome reference build name | "GRCh38" |
| File | ref_dict | Yes | The genome reference dict file | |
| File | ref_fasta | Yes | The genome reference fasta file | |
| File | ref_fasta_index | Yes | The genome reference fasta index file | |
| File | scattered_calling_intervals_list | No | File containing list of scattered calling interval files for Haplotype Caller | |
| File | cov_roi_bed | No | The BED that defines the region of interest for coverage analysis; see DepthOfCoverage.md for details | |
| Array[RoiAndRefGeneFilePair] | cov_roi_genes | No | List of ROI and ref gene file pairs; see DepthOfCoverage.md for details | |
| File | cov_gene_names | No | Tab-delimited file of gene information; see DepthOfCoverage.md for details | |
| String | cov_docker_image | Yes | The name of the coverage Docker image to run the coverage summary task | |
| String | gatk3_docker_image | No | The name of the GATK3 Docker image for coverage analysis | "broadinstitute/gatk3:3.7-0" |
| File | target_roi_bed | Yes | THe BED that defines the target region of interest for annotation and filtration | |
| File | alamut_db | Yes | The database file for Alamut batch | |
| File | alamut_fields_tsv | No | The file that defines how the Alamut output is transformed back to a VCF | |
| String | alamut_db_name | No | The database name for the Alamut batch ini file | "alamut_db" |
| String | alamut_server | No | The server name for the Alamut batch ini file | "a-ht-na.interactive-biosoftware.com" |
| String | alamut_port | No | The server port for the Alamut batch ini file | "80" |
| String | alamut_user_secret_name | No | The GCP secret name that contains the user stanza for the Alamut batch ini file | "alamut-batch-ini-user" |
| Int | alamut_queue_limit | No | The maximum number of concurrent Alamut batch processes permitted | 4 |
| String | alamut_queue_folder | No | The shared storage location for Alamut concurrency management | "gs://biofx-task-queue/alamut" |
| String | alamut_docker_image | Yes | The name of the Alamut Docker image using to run Alamut Batch task | |
| Boolean | alamut_save_working_files | No | Whether or not to retain intermediate Alamut Batch task files | false |
| String | alamut_anno_src_id | No | When removing already annotated variants prior to Alamut annotation, the annotation source id to query for | "228" |
| String | alamut_anno_min_age | No | When removing already annotated variants prior to Alamut annotation, the annotation minimum timestamp to query for (ISO8601 duration) | "P6M" |
| String | qceval_project_type | No | The type of rules to apply for the QC evaluation task, one of "WGS", "WGS_DRAGEN", "WES" or "NONE" | "WGS" |
| String | qceval_docker_image | Yes | The name of the Docker image to run the QC evaluation task | |
| File | gnomad_coverage_file | No | The gnomad coverage data file | |
| File | gnomad_coverage_file_idx | No | The gnomad coverage data index file | |
| Array[String] | gnomad_headers | No | List of VCF headers to add when annotating VCF with gnomad coverage data.  For example, ##INFO=<ID=DP_gnomadG,Number=1,Type=Float,Description="Read depth of GnomAD Genome"> | |
| String | gnomad_column_list | No | The column list to pass to bcftools annotate for gnomad coverage annotation | "CHROM,POS,INFO/DP_gnomadG" |
| String | sample_data_load_config_name | Yes | The FAST load configuration name for the sample data VCF | |
| String | gnomad_data_load_config_name | Yes | The FAST load configuration name for the gnomad coverage VCF | |
| String | alamut_data_load_config_name | Yes | The FAST load configuration name for the Alamut annotated VCF | |
| Array[String] | fast_annotated_sample_data_regions | No | The list of regions to include in the FAST annotated sample data; each element is a "name:applyMask" pair | |
| Array[String] | fast_annotated_sample_data_scripts | No | The list of custom scripts to run on the FAST annotated sample data after creation | |
| String | fast_annotated_sample_data_saved_filter_name | No | The saved filter to apply to the FAST annotated sample data | |
| String | igvreport_docker_image | Yes | The name of the Docker image to run the IGV report task | |
| String | fast_parser_image | Yes | The name of the Docker image to run the FAST output parser task | |
| File | gil_transcript_exon_count | Yes | A tab delimited file of transcript id and exon count |
| String | fast_parser_sample_type | Yes | The sample type flag for the FAST output parser: S for single-sample Exome or M for multi-sample Exome or B for batch/Biobank or N for NVA-Lite |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | cov_wgs_sample_summary | Always | Coverage metrics summary for all bases |
| File | cov_wgs_sample_statistics | Always |  Coverage metrics statistics for all bases |
| File | cov_roi_sample_interval_summary | If coverage ROI is specified | Coverage metrics interval summary for region of interest |
| File | cov_roi_sample_interval_statistics | If coverage ROI is specified |  Coverage metrics interval statistics for region of interest |
| File | cov_roi_sample_statistics | If coverage ROI is specified | Coverage metrics statistics for region of interest |
| File | cov_roi_sample_summary | If coverage ROI is specified | Coverage metrics summary for region of interest |
| File | cov_roi_sample_cumulative_coverage_counts | If coverage ROI is specified | Coverage cumulative counts for the region of interest |
| File | cov_roi_sample_cumulative_coverage_proportions | If coverage ROI is specified | Coverage cumulative proportions for the region of interest |
| File | cov_mt_summary | If coverage ROI genes are specified | Mitochondrial gene coverage summary |
| File | cov_gene_summary | If coverage ROI genes are specified | Gene coverage summary |
| File | cov_gene_summary_unknown | If coverage ROI genes are specified | Unknown entries from the gene coverage summary |
| File | cov_gene_summary_entrez | If coverage ROI genes are specified | Gene coverage summary enriched with Entrez IDs |
| File | vcf | If variant calling is run | Variants called from BAM/CRAM |
| File | target_vcf_gz | Always | VCF file filtered to the target region of interest |
| File | alamut_vcf_gz | Always | Target VCF file with Alamut annotations |
| File | qceval_vcf_gz | Always | Target VCF file annotated with QC Evaluation |
| File | gnomad_vcf_gz | If gnomad coverage data is specified | Target VCF annotated with gnomad coverage data |
| File | fast_export_file | Always | Tab-delimited export of annotated sample data from FAST |
| File | fast_summary_file | Always | Summary of FAST processing parameters |
| File | igv_report | Always | HTML-based IGV report file |
| File | fast_parsed_output | Always | Parsed FAST export |
| File | nva_report | Always | NVA report Excel document |