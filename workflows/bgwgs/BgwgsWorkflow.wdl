version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"
import "../../steps/DepthOfCoverage.wdl"
import "../../steps/AlamutBatch.wdl"
import "../../steps/VCFUtils.wdl"
import "../../steps/FASTUtils.wdl"
import "../../steps/QCEval.wdl"
import "../../steps/IgvReport.wdl"
import "../../steps/FASTOutputParser.wdl"
import "../../steps/Utilities.wdl"

workflow BgwgsWorkflow {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image
        # bcftools docker image
        String bcftools_docker_image
        # subject, sample id and data location
        String subject_id
        String sample_id
        String sample_data_location
        Boolean fetch_cram = true
        Boolean fetch_bam = true
        Boolean do_variant_calling = true
        Int fetch_disk_size = 75
        # reference genome files
        String reference_build = "GRCh38"
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        # variant calling inputs
        File? scattered_calling_intervals_list
        # depth of coverage inputs
        File? cov_roi_bed
        Array[RoiAndRefGeneFilePair?] cov_roi_genes
        File? cov_gene_names
        String cov_docker_image
        String gatk3_docker_image = "broadinstitute/gatk3:3.7-0"
        # vcf filter inputs
        File target_roi_bed
        # alamut inputs
        File alamut_db
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        String alamut_docker_image
        Boolean alamut_save_working_files = false
        String alamut_anno_src_id = "228"
        String alamut_anno_min_age = "P6M"
        # qceval inputs
        String qceval_project_type = "WGS"
        String qceval_docker_image
        # gnomad annotation inputs
        File? gnomad_coverage_file
        File? gnomad_coverage_file_idx
        Array[String]? gnomad_headers
        String gnomad_column_list = "CHROM,POS,INFO/DP_gnomadG"
        # FAST loading inputs
        Boolean has_haploid_sites
        String sample_data_load_config_name
        String gnomad_data_load_config_name
        String alamut_data_load_config_name
        Array[String]? fast_annotated_sample_data_regions
        Array[String]? fast_annotated_sample_data_scripts
        String? fast_annotated_sample_data_saved_filter_name
        # Reporting steps
        String igvreport_docker_image
        String fast_parser_image
        File gil_transcript_exon_count
        String fast_parser_sample_type
    }

    # Prefer a BAM file to avoid conversion
    if (fetch_bam) {
        call FileUtils.FetchFilesTask as FetchBam {
            input:
                data_location = sample_data_location,
                file_types = ["bam"],
                recursive = false,
                file_match_keys = [subject_id, sample_id],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
        }
    }

    # no BAM, fallback to CRAM
    if (fetch_cram && (!defined(FetchBam.bam) || !defined(FetchBam.bai))) {
        call FileUtils.FetchFilesTask as FetchCram {
            input:
                data_location = sample_data_location,
                file_types = ["cram"],
                recursive = false,
                file_match_keys = [subject_id, sample_id],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
        }
    }

    # Validate that we got workable files from the data location
    if (!defined(FetchBam.bam) && !defined(FetchCram.cram)) {
        call Utilities.FailTask as MissingBamOrCramFailure {
            input:
                error_message = "BAM or CRAM file not found in ~{sample_data_location}"
        }
    }
    if (defined(FetchCram.cram) && !defined(FetchCram.crai)) {
        call Utilities.FailTask as MissingCraiFailure {
            input:
                error_message = "Index file for CRAM " + basename(select_first([FetchCram.cram])) + " not found"
        }
    }
    if (!defined(FetchCram.cram) && defined(FetchBam.bam) && !defined(FetchBam.bai)) {
        call Utilities.FailTask as MissingBaiFailure {
            input:
                error_message = "Index file for BAM " + basename(select_first([FetchCram.bam])) + " not found"
        }
    }

    # Convert CRAM to BAM
    if (defined(FetchCram.cram) && defined(FetchCram.crai)) {
        call HaplotypeCallerGvcfGATK4.CramToBamTask {
            input:
                input_cram = select_first([FetchCram.cram]),
                sample_name = basename(select_first([FetchCram.cram]), ".cram"),
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
                samtools_path = "samtools"
        }
    }

    # final inputs - either CRAM or BAM
    File sample_bam = select_first([FetchBam.bam, CramToBamTask.output_bam])
    File sample_bai = select_first([FetchBam.bai, CramToBamTask.output_bai])

    # Run haplotype caller
    if (do_variant_calling) {
        call HaplotypeCallerGvcfGATK4.HaplotypeCallerGvcf_GATK4 {
            input:
                input_bam = sample_bam,
                input_bam_index = sample_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                scattered_calling_intervals_list = select_first([scattered_calling_intervals_list]),
                make_gvcf = false
        }
    }
    if (!do_variant_calling) {
        # Fetch the VCF files from the storage location
        call FileUtils.FetchFilesTask as FetchVcf {
            input:
                data_location = sample_data_location,
                file_types = ["vcf"],
                recursive = false,
                file_match_keys = [subject_id, sample_id],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name
        }
        if (!defined(FetchVcf.vcf)) {
            call Utilities.FailTask as MissingVcfFailure {
                input:
                    error_message = "VCF file not found in ~{sample_data_location}"
            }
        }
    }

    File sample_vcf = select_first([HaplotypeCallerGvcf_GATK4.output_vcf, FetchVcf.vcf])

    # Run depth of coverage
    call DepthOfCoverage.DepthOfCoverageWorkflow {
        input:
            run_wgs = true,
            output_basename = subject_id + "_" + sample_id + ".cov",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            bam = sample_bam,
            bai = sample_bai,
            roi_all_bed = cov_roi_bed,
            roi_genes = cov_roi_genes,
            gene_names = cov_gene_names,
            cov_docker_image = cov_docker_image,
            gatk_docker_image = gatk3_docker_image
    }

    # Filter called VCF to target region
    call VCFUtils.FilterVCFWithBEDTask {
        input:
            input_vcf = sample_vcf,
            input_bed = target_roi_bed,
            output_basename = subject_id + "_" + sample_id + ".target",
            docker_image = "~{bcftools_docker_image}"
    }

    # Annotate target VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
            project_type = qceval_project_type,
            output_basename = subject_id + "_" + sample_id + ".qceval",
            docker_image = qceval_docker_image
    }

    # Pre-filter the VCF file to input to Alamut to remove
    #  variants that already have Alamut annotations
    call FASTUtils.FASTRemoveAlreadyAnnotatedFromVCFTask {
        input:
            input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
            output_basename = subject_id + "_" + sample_id + ".alamutprefilter",
            reference_build = reference_build,
            annotation_source_id = alamut_anno_src_id,
            annotation_min_age = alamut_anno_min_age,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    # Annotate target VCF with Alamut
    call AlamutBatch.AlamutBatchTask {
        input:
            input_vcf = FASTRemoveAlreadyAnnotatedFromVCFTask.output_vcf_gz,
            output_basename = subject_id + "_" + sample_id + ".alamut",
            alamut_db = alamut_db,
            alamut_fields_tsv = alamut_fields_tsv,
            alamut_db_name = alamut_db_name,
            alamut_server = alamut_server,
            alamut_port = alamut_port,
            reference_build = reference_build,
            gcp_project_id = gcp_project_id,
            alamut_user_secret_name = alamut_user_secret_name,
            alamut_queue_limit = alamut_queue_limit,
            alamut_queue_folder = alamut_queue_folder,
            docker_image = alamut_docker_image,
            output_working_files = alamut_save_working_files
    }

    # Annotated target VCF with Gnomad coverage
    if (defined(gnomad_coverage_file) && defined(gnomad_headers) && defined(gnomad_coverage_file_idx)) {
        call VCFUtils.AnnotateVCFTask as AnnotateGnomadTask {
            input:
                input_vcf = FilterVCFWithBEDTask.output_vcf_gz,
                output_basename = subject_id + "_" + sample_id + ".gnomad",
                annotations_file = select_first([gnomad_coverage_file]),
                annotations_idx_file = select_first([gnomad_coverage_file_idx]),
                headers_file = write_lines(select_first([gnomad_headers])),
                column_list = gnomad_column_list,
                docker_image = "~{bcftools_docker_image}"
        }
    }

    # Load sample data to FAST
    call FASTUtils.FASTDataLoadTask as QCEvalLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = QCEvalTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            sample_data_name = subject_id + "_" + sample_id,
            lab_batch_name = sample_id,
            data_load_config_name = sample_data_load_config_name,
            data_load_target = "SAMPLE_DATA",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    # Load Alamut annotations to FAST
    call FASTUtils.FASTDataLoadTask as AlamutLoadTask {
        input:
            reference_build = reference_build,
            vcf_file = AlamutBatchTask.output_vcf_gz,
            has_haploid_sites = has_haploid_sites,
            data_load_config_name = alamut_data_load_config_name,
            data_load_target = "ANNOTATION_DATA",
            merge_strategy = "MERGE",
            annotation_record_ts = "now",
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    # Load Gnomad annotations to FAST
    if (defined(AnnotateGnomadTask.output_vcf_gz)) {
        call FASTUtils.FASTDataLoadTask as GnomadLoadTask {
            input:
                reference_build = reference_build,
                vcf_file = select_first([AnnotateGnomadTask.output_vcf_gz]),
                has_haploid_sites = has_haploid_sites,
                data_load_config_name = gnomad_data_load_config_name,
                data_load_target = "ANNOTATION_DATA",
                merge_strategy = "MERGE",
                annotation_record_ts = "now",
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
    }

    # Wait for data loads to complete
    call FASTUtils.FASTWaitForDataLoadsTask {
        input:
            data_load_ids = select_all([QCEvalLoadTask.data_load_id, AlamutLoadTask.data_load_id, GnomadLoadTask.data_load_id]),
            check_interval_secs = 300,
            max_checks = 72,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.completed) {
        call Utilities.FailTask as DataLoadTimeout {
            input:
                error_message = "One or more FAST data loads did not complete within allotted time"
        }
    }
    if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.succeeded) {
        call Utilities.FailTask as DataLoadFailure {
            input:
                error_message = "One or more FAST data loads did not succeed"
        }
    }

    # Create annotated sample data
    if (FASTWaitForDataLoadsTask.wait_result.total == FASTWaitForDataLoadsTask.wait_result.succeeded) {
        call FASTUtils.FASTCreateAnnotatedSampleDataTask {
            input:
                annotated_sample_data_name = subject_id + "_" + sample_id,
                sample_data_names_and_labels = [subject_id + "_" + sample_id],
                region_names_and_masks = fast_annotated_sample_data_regions,
                scripts = fast_annotated_sample_data_scripts,
                saved_filter_name = fast_annotated_sample_data_saved_filter_name,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }

        # Wait for annotated sample data to complete
        call FASTUtils.FASTWaitForAnnotatedSampleDataTask {
            input:
                annotated_sample_data_name = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output,
                check_interval_secs = 300,
                max_checks = 72,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
        if (FASTWaitForAnnotatedSampleDataTask.wait_result.total != FASTWaitForAnnotatedSampleDataTask.wait_result.completed) {
            call Utilities.FailTask as AnnotatedSampleTimeout {
                input:
                    error_message = "FAST annotated sample data creation did not complete within allotted time"
            }
        }
        if (FASTWaitForAnnotatedSampleDataTask.wait_result.total != FASTWaitForAnnotatedSampleDataTask.wait_result.succeeded) {
            call Utilities.FailTask as AnnotatedSampleFailure {
                input:
                    error_message = "FAST annotated sample data creation did not succeed"
            }
        }

        if (FASTWaitForAnnotatedSampleDataTask.wait_result.total == FASTWaitForAnnotatedSampleDataTask.wait_result.succeeded) {
            call FASTUtils.FASTExportAnnotatedSampleDataTask {
                input:
                    annotated_sample_data_name = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output,
                    format = "TXT",
                    output_basename = subject_id + "_" + sample_id + ".fastexport",
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    docker_image = orchutils_docker_image
            }
            call FASTOutputParser.FASTOutputParserTask {
                input:
                    fast_output_file = FASTExportAnnotatedSampleDataTask.output_file,
                    sample_type = fast_parser_sample_type,
                    reference_build = reference_build,
                    oms_query = "Y",
                    transcript_exonNum = gil_transcript_exon_count,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    fast_parser_image = fast_parser_image
            }
            if (defined(FASTOutputParserTask.parsed_report)) {
                call IgvReport.IgvReportFromParsedFASTOutputTask {
                    input:
                        bam_cram = sample_bam,
                        bai_crai = sample_bai,
                        parsed_fast_output = select_first([FASTOutputParserTask.parsed_report]),
                        output_basename = subject_id + "_" + sample_id + ".igvreport",
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        docker_image = igvreport_docker_image
                }
            }
            call BgwgsFASTSummaryTask {
                input:
                    subject_id = subject_id,
                    sample_id = sample_id,
                    fast_annotated_sample_data_regions = fast_annotated_sample_data_regions,
                    fast_annotated_sample_data_saved_filter_name = fast_annotated_sample_data_saved_filter_name
            }
        }
    }

    output {
        # coverage outputs
        File? cov_wgs_sample_summary = DepthOfCoverageWorkflow.wgs_sample_summary
        File? cov_wgs_sample_statistics = DepthOfCoverageWorkflow.wgs_sample_statistics
        File? cov_roi_sample_interval_summary = DepthOfCoverageWorkflow.roi_sample_interval_summary
        File? cov_roi_sample_interval_statistics = DepthOfCoverageWorkflow.roi_sample_interval_statistics
        File? cov_roi_sample_statistics = DepthOfCoverageWorkflow.roi_sample_statistics
        File? cov_roi_sample_summary = DepthOfCoverageWorkflow.roi_sample_summary
        File? cov_roi_sample_cumulative_coverage_counts = DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_counts
        File? cov_roi_sample_cumulative_coverage_proportions = DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_proportions
        File? cov_mt_summary = DepthOfCoverageWorkflow.mt_summary
        File? cov_gene_summary = DepthOfCoverageWorkflow.gene_summary
        File? cov_gene_summary_unknown = DepthOfCoverageWorkflow.gene_summary_unknown
        File? cov_gene_summary_entrez = DepthOfCoverageWorkflow.gene_summary_entrez
        # haplotype caller output
        File? vcf = HaplotypeCallerGvcf_GATK4.output_vcf
        # filtered VCF
        File target_vcf_gz = FilterVCFWithBEDTask.output_vcf_gz
        # annotated VCFs
        File alamut_vcf_gz = AlamutBatchTask.output_vcf_gz
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        File? gnomad_vcf_gz = AnnotateGnomadTask.output_vcf_gz
        # FAST export file
        File? fast_export_file = FASTExportAnnotatedSampleDataTask.output_file
        # FAST summary file
        File? fast_summary_file = BgwgsFASTSummaryTask.fast_summary_file
        # IGV report
        File? igv_report = IgvReportFromParsedFASTOutputTask.igv_report
        # FAST Parsed output and NVA report
        File? fast_parsed_output = FASTOutputParserTask.parsed_report
        File? nva_report = FASTOutputParserTask.nva_report
    }
}

task BgwgsFASTSummaryTask {
    input {
        String subject_id
        String sample_id
        Array[String]? fast_annotated_sample_data_regions
        String? fast_annotated_sample_data_saved_filter_name
    }

    command <<<
        summary_file="~{subject_id}_~{sample_id}_FAST_summary.txt"

        printf "PM_number\tSaved_filter\tRegions\n" >> "${summary_file}"
        printf "%s\t%s\t%s\n" "~{subject_id}" "~{fast_annotated_sample_data_saved_filter_name}" "~{sep=',' fast_annotated_sample_data_regions}" >> "${summary_file}"
    >>>

    runtime {
        docker: "ubuntu:22.04"
    }

    output {
        File fast_summary_file = subject_id + "_" + sample_id + "_FAST_summary.txt"
    }
}