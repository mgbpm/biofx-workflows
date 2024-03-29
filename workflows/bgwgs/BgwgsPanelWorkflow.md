# BGWGS Panel Workflow
The BGWGS (bigwigs) Genome Panels Workflow starts with a single CRAM or BAM file and provides variant calling, and LMM reporting capabilities based on a specific panel.

## Input Parameters
 Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | gcp_project_id | No | The GCP project to fetch secrets from | "mgb-lmm-gcp-infrast-1651079146" |
| String | workspace_name | Yes | The name of the current workspace (for secret retrieval) | |
| String | orchutils_docker_image | No | The name of the orchestration utils Docker image for FAST and file movement tasks | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230828" |
| String | genome_panels_docker_image | No | The name of the genome panels Docker image for VCF annotation | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/genome-panels:20240208" |
| String | subject_id | Yes | The sample id associated with the data (PM-ID) | |
| String | sample_id | Yes | The sample Run id associated with the data (SM-ID) | |
| String | sample_lmm_id | Yes | The sample LMM id associated with the data (LMM-ID) | |
| String | sample_barcode | Yes | The sample barcode associated with the data | |
| String | batch_id | Yes | The batch id associated with the data | |
| String | sample_data_location | Yes | The cloud storage URL where the sample source data is located | |
| Boolean | fetch_cram | No | Whether or not to fetch the CRAM (primarily for testing) | true |
| Array[String] | fetch_cram_filter_keys | No | The list of strings that must appear in the CRAM file path | [subject_id, sample_id] |
| Array[FileMatcher]? | fetch_cram_file_matchers | No | The list of file matchers to use for CRAM file fetching | |
| Boolean | fetch_bam | No | Whether or not to fetch the BAM (primarily for testing) | true |
| Array[String] | fetch_bam_filter_keys | No | The list of strings that must appear in the BAM file path | [subject_id, sample_id] |
| Array[FileMatcher]? | fetch_bam_file_matchers | No | The list of file matching rules for BAM fetching | |
| Array[String] | fetch_vcf_filter_keys | No | The list of strings that must appear in the VCF file path | [subject_id, sample_id] |
| Array[FileMatcher]? | fetch_vcf_file_matchers | No | The list of file matching rules for VCF fetching | |
| Boolean | fetch_files_verbose | No | If true, generate verbose output from file fetch tasks | false |
| String | test_code | Yes | The testcode or panel to be run | |
| File | ref_dict | Yes | The genome reference dict file | |
| File | ref_fasta | Yes | The genome reference fasta file | |
| File | ref_fasta_index | Yes | The genome reference fasta index file | |
| File | bait_file | Yes | File containing list of bait intervals | |
| File | extra_amplicons_file | Yes | File containing list of extra amplicons | |
| File | duplicate_amplicons_file | Yes | File containing list of duplicate amplicons | |
| Boolean | do_variant_calling | No | Whether or not to do variant calling | true |
| File | target_intervals | Yes | File containing list of interval files for Haplotype Caller | |
| File | dbsnp | No | HG38 DBSNP vcf file | |
| File | dbsnp_vcf_index | No | HG38 DBSNP vcf index file | |
| String | igvreport_docker_image | No | The name of the IGV Docker image for the IGV report generation task | "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511" |


## Output Parameters
| Type | Name | Description |
| :--- | :--- | :--- |
| File | all_calls_vcf_file | VCF file with all calls |
| File | all_calls_vcf_idx_file | Index file for VCF file with all calls |
| File | gvcf_file | GVCF file from the Variant Calling task |
| File | gvcf_idx_file | Index file for GVCF file from the Variant Calling task |
| File | ref_positions_vcf_file | Ref positions file for sorting the all calls file |
| File | all_bases_vcf_file | VCF file with all bases |
| File | all_bases_idx_vcf_file | Index file for VCF file with all bases |
| File | all_bases_noChr_vcf_file | VCF file with all bases without the CHR prefix |
| File | snps_out_file | Text file with common YSNPs |
| File | xls_report_out_file | LMM Variant Detection report in XLS format |
| File | xml_report_out_file | LMM Variant Detection report in XML format |
| File | follow_up_regions_file | Follow up regions BED file used to generate IGV report |
| File | lmm_igv_report | LMM IGV report from the follow up BED file |
