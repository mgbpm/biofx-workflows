version 1.0

import "../../steps/Utilities.wdl" # for testing inputs
import "../../steps/FileUtils.wdl" # for fetching files
import "../../steps/VCFUtils.wdl" # for processing input data
import "../../steps/AlamutBatch.wdl" # for Alamut annotation
import "../../steps/QCEval.wdl" # for qc of each vcf
import "../../steps/FASTUtils.wdl" # for loading to FAST
import "../../steps/FASTOutputParser.wdl" # for output
import "../pgxrisk/PGxWorkflow.wdl" # for running PGx workflow
import "../pgxrisk/RiskAllelesWorkflow.wdl" # for running Risk workflow

workflow BahrainPipelinesWorkflow {
    input {
        # Data prep inputs
        Array[String] sample_ids
        Array[String] sample_data_locations
        String batch_name
        File target_roi_bed
        String pipeline_to_run # either "monogenic" or "screening"
        # Python docker image
        String python_docker_image = "python:3.10"
        # Bcftools docker image
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker image
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129"
        # Reference genome files
        String reference_build = "GRCh38"
        File ref_fasta
        File ref_fasta_index
        File? ref_dict
        File? dbsnp_vcf
        File? dbsnp_vcf_index
        # PGx inputs
        String pgx_test_code = "lmPGX-pnlC_L"
        String pgx_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20230711"
        File pgx_workflow_fileset = "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_files-20220118.tar"
        File pgx_roi_bed = "gs://lmm-reference-data/pgx/lmPGX-pnlC_L_genotyping-chr-20220118.bed"
        # Risk alleles inputs
        String risk_alleles_test_code = "lmRISK-pnlB_L"
        String risk_alleles_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20230724"
        File risk_alleles_workflow_fileset = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar"
        File risk_alleles_roi_bed = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed"
        # Alamut inputs
        File alamut_db
        File? alamut_fields_tsv
        String alamut_db_name = "alamut_db"
        String alamut_server = "a-ht-na.interactive-biosoftware.com"
        String alamut_port = "80"
        String alamut_user_secret_name = "alamut-batch-ini-user"
        Int alamut_queue_limit = 4
        String alamut_queue_folder = "gs://biofx-task-queue/alamut"
        String alamut_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/alamut:20230630"
        Boolean alamut_save_working_files = false
        String alamut_anno_src_id = "228"
        String alamut_anno_min_age = "P6M"
        # QCEval inputs
        String qceval_project_type = "WGS"
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005"
        # FAST inputs
        Boolean has_haploid_sites = false
        String sample_data_load_config_name = "Sample_VCF_PPM_Eval"
        String alamut_data_load_config_name = "Alamut"
        Array[String]? fast_annotated_sample_data_regions
        Array[String]? fast_annotated_sample_data_scripts
        String fast_annotated_sample_data_saved_filter_name
        Int fast_data_load_wait_interval_secs = 300
        Int fast_data_load_wait_max_intervals = 72
        Int fast_adi_wait_interval_secs = 600
        Int fast_adi_wait_max_intervals = 144
        String fast_parser_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/fastoutputparser:20231206"
        File gil_transcript_exon_count = "gs://lmm-reference-data/annotation/gil_lmm/transcript_exonNum.txt"
        String fast_parser_sample_type
    }

    ## Test that inputs match the pipeline being run
    # Test pipeline to run -- should be monogenic or screening
    if ((pipeline_to_run != "monogenic") && (pipeline_to_run != "screening")) {
        String PipelineNameFail = "Pipeline to run should be either 'monogenic' or 'screening.'"
    }
    # Test that the number of sample IDs matches the number of sample data locations
    if (length(sample_ids) != length(sample_data_locations)) {
        String SampleArraysLengthFail = "The number of sample IDs does not match the number of sample data locations."
    }
    if (pipeline_to_run == "monogenic") {
        # Test that there is not a large number of samples run for the monogenic pipeline
        if (length(sample_ids) > 50) {
            String MonogenicFileNumberFail = "To run the monogenic pipeline, the number of input files must be less than 50."
        }
        # Test that the FAST parser sample type matches the pipeline being run
        if (fast_parser_sample_type != "M") {
            String MonogenicParserSampleTypeFail = "When running the monogenic pipeline, the FAST parser sample type should be 'M'."
        }
    }
    if (pipeline_to_run == "screening") {
        # Test that the FAST parser sample type matches the pipeline being run
        if (fast_parser_sample_type != "B") {
            String ScreeningParserSampleTypeFail = "When running the screening pipeline, the FAST parser sample type should be 'B'."
        }
        # Test that there are all inputs for the PGx/Risk workflows
        if (!defined(ref_dict) || !defined(dbsnp_vcf) || !defined(dbsnp_vcf_index)) {
            String ScreeningRefInputFail = "When running the screening pipeline, include inputs for ref_dict, dbsnp_vcf, and dbsnp_vcf_index."
        }
    }
    String input_error_message = select_first([PipelineNameFail, MonogenicFileNumberFail, SampleArraysLengthFail, MonogenicParserSampleTypeFail, ScreeningParserSampleTypeFail, ScreeningRefInputFail, ""])
    if (input_error_message != "") {
        call Utilities.FailTask as InputParameterError {
            input:
                error_message = input_error_message
        }
    }
    
    ## Fetch all files (CRAM, VCF, and index files)
    if (input_error_message == "") {
        scatter (i in range(length(sample_ids))) {
            # Fetch sample VCFs and CRAMs for the screening pipeline
            if (pipeline_to_run == "screening") {
                call FileUtils.FetchFilesTask as FetchScreeningFiles {
                    input:
                        data_location = sample_data_locations[i],
                        recursive = true,
                        file_types = [ "vcf", "cram" ],
                        file_match_keys = [ sample_ids[i] ],
                        docker_image = orchutils_docker_image,
                        disk_size = 50,
                        gcp_project_id = gcp_project_id,
                        workspace_name = workspace_name
                }
                if (!defined(FetchScreeningFiles.cram)) {
                    call Utilities.FailTask as CRAMNotFound {
                        input:
                            error_message = "CRAM for sample " + sample_ids[i] + " not found in " + sample_data_locations[i]
                    }
                }
                if (!defined(FetchScreeningFiles.crai)) {
                    call Utilities.FailTask as CRAMIndexNotFound {
                        input:
                            error_message = "CRAM index for sample " + sample_ids[i] + " not found in " + sample_data_locations[i]
                    }
                }
            }
            # Fetch only VCFs for the monogenic pipeline
            if (pipeline_to_run == "monogenic") {
                call FileUtils.FetchFilesTask as FetchMonogenicFiles {
                    input:
                        data_location = sample_data_locations[i],
                        recursive = true,
                        file_types = [ "vcf" ],
                        file_match_keys = [ sample_ids[i] ],
                        docker_image = orchutils_docker_image,
                        disk_size = 20,
                        gcp_project_id = gcp_project_id,
                        workspace_name = workspace_name
                }
            }
            if (!defined(select_first([FetchScreeningFiles.vcf, FetchMonogenicFiles.vcf]))) {
                call Utilities.FailTask as VCFFileNotFound {
                    input:
                        error_message = "VCF file for sample " + sample_ids[i] + " not found in " + sample_data_locations[i]
                }
            }
            if (!defined(select_first([FetchScreeningFiles.vcf_index, FetchMonogenicFiles.vcf_index]))) {
                call Utilities.FailTask as VCFIndexNotFound {
                    input:
                        error_message = "VCF index for sample " + sample_ids[i] + " not found in " + sample_data_locations[i]
                }
            }
            File found_vcf = select_first([FetchScreeningFiles.vcf, FetchMonogenicFiles.vcf])
            File found_vcf_idx = select_first([FetchScreeningFiles.vcf_index, FetchMonogenicFiles.vcf_index])
        }
    }
    # Coerce object types to Array[File] for future tasks
    Array[File] fetched_vcfs = select_first([found_vcf])
    Array[File] fetched_vcfs_idx = select_first([found_vcf_idx])
    Array[File] fetched_crams = select_all(select_first([FetchScreeningFiles.cram]))
    Array[File] fetched_crams_idx = select_all(select_first([FetchScreeningFiles.crai]))

    ## Run PGx & Risk workflows with fetched CRAMs
    if (pipeline_to_run == "screening") {
        scatter (i in range(length(fetched_crams))) {
            call PGxWorkflow.PGxWorkflow as PGxWorkflowAlias {
                input:
                    input_cram = fetched_crams[i],
                    input_crai = fetched_crams_idx[i],
                    sample_id = sample_ids[i],
                    accession_id = sample_ids[i],
                    test_code = pgx_test_code,
                    reference_fasta = ref_fasta,
                    reference_fasta_fai = ref_fasta_index,
                    reference_dict = select_first([ref_dict]),
                    roi_bed = pgx_roi_bed,
                    dbsnp = select_first([dbsnp_vcf]),
                    dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                    workflow_fileset = pgx_workflow_fileset,
                    mgbpmbiofx_docker_image = pgx_docker_image
            }
            call RiskAllelesWorkflow.RiskAllelesWorkflow as RiskAllelesWorkflowAlias {
                input:
                    input_cram = fetched_crams[i],
                    input_crai = fetched_crams_idx[i],
                    sample_id = sample_ids[i],
                    accession_id = sample_ids[i],
                    test_code = risk_alleles_test_code,
                    reference_fasta = ref_fasta,
                    reference_fasta_fai = ref_fasta_index,
                    reference_dict = select_first([ref_dict]),
                    roi_bed = risk_alleles_roi_bed,
                    dbsnp = select_first([dbsnp_vcf]),
                    dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                    workflow_fileset = risk_alleles_workflow_fileset,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    mgbpmbiofx_docker_image = risk_alleles_docker_image
            } 
        }
    }

    ## Prep data -- filter vcfs, merge them, and create collective vcf (if necessary)
    # If a target ROI bed file is given, filter the single vcfs to the bed file
    scatter (i in range(length(fetched_vcfs))) {
        call PrepSampleVCFTask {
            input:
                input_vcf = fetched_vcfs[i],
                input_vcf_idx = fetched_vcfs_idx[i],
                target_roi_bed = target_roi_bed,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                docker_image = bcftools_docker_image
        }
    }
    # Create merged VCF from filtered sample files
    call VCFUtils.MergeVCFsTask as MergedVCF {
        input:
            input_vcfs = PrepSampleVCFTask.output_vcf_gz,
            input_vcfs_idx = PrepSampleVCFTask.output_vcf_idx,
            sorted = true,
            output_basename = batch_name + ".merged",
            docker_image = bcftools_docker_image
    }
    # If running the screening pipeline, create a collective VCF
    if (pipeline_to_run == "screening") {
        call VCFUtils.MakeCollectiveVCFTask as CollectiveVCF {
            input:
                input_vcf = MergedVCF.output_vcf_gz,
                docker_image = python_docker_image
        }
    }

    ## Annotate with Alamut and load Alamut annotations to FAST
    # Pre-filter the collective/merged VCF to input to Alamut to remove
    #  variants that already have Alamut annotations
    call FASTUtils.FASTRemoveAlreadyAnnotatedFromVCFTask {
        input:
            input_vcf = select_first([CollectiveVCF.output_vcf_gz, MergedVCF.output_vcf_gz]),
            output_basename = batch_name + "_all.alamutprefilter",
            reference_build = reference_build,
            annotation_source_id = alamut_anno_src_id,
            annotation_min_age = alamut_anno_min_age,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    # Annotate collective or merged VCF with Alamut
    call AlamutBatch.AlamutBatchTask {
        input:
            input_vcf = FASTRemoveAlreadyAnnotatedFromVCFTask.output_vcf_gz,
            output_basename = batch_name + "_all.alamut",
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
    call FASTUtils.FASTWaitForDataLoadsTask as WaitForAlamutLoad {
        input:
            data_load_ids = [AlamutLoadTask.data_load_id],
            check_interval_secs = fast_data_load_wait_interval_secs,
            max_checks = fast_data_load_wait_max_intervals,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }

    ## Annotate each individual sample VCF with QCEval INFO field and load annotations into FAST
    scatter (file in PrepSampleVCFTask.output_vcf_gz) {
        String sample_name = sub(sub(basename(file), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", ""), ".filtered.norm", "")
        String qc_fast_sample_data_name = sample_name + "_" + batch_name
        # Annotate individual VCF with QCEval INFO field
        call QCEval.QCEvalTask {
            input:
                input_vcf = file,
                project_type = qceval_project_type,
                output_basename = sample_name + ".qceval",
                docker_image = qceval_docker_image
        }
        # Load sample data to FAST
        call FASTUtils.FASTDataLoadTask as QCEvalLoadTask {
            input:
                reference_build = reference_build,
                vcf_file = QCEvalTask.output_vcf_gz,
                has_haploid_sites = has_haploid_sites,
                sample_data_name = qc_fast_sample_data_name,
                data_load_config_name = sample_data_load_config_name,
                data_load_target = "SAMPLE_DATA",
                annotation_record_ts = "now",
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
    }

    ## Wait for both data loads to complete
    call FASTUtils.FASTWaitForDataLoadsTask as WaitForQCLoads {
        input:
            data_load_ids = QCEvalLoadTask.data_load_id,
            check_interval_secs = 300,
            max_checks = 72,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = orchutils_docker_image
    }
    if ((WaitForQCLoads.wait_result.total != WaitForQCLoads.wait_result.completed) || (WaitForAlamutLoad.wait_result.total != WaitForAlamutLoad.wait_result.completed)) {
        call Utilities.FailTask as DataLoadTimeout {
            input:
                error_message = "One or more FAST data loads did not complete within allotted time"
        }
    }
    if ((WaitForQCLoads.wait_result.total != WaitForQCLoads.wait_result.succeeded) || (WaitForAlamutLoad.wait_result.total != WaitForAlamutLoad.wait_result.succeeded)) {
        call Utilities.FailTask as DataLoadFailure {
            input:
                error_message = "One or more FAST data loads did not succeed"
        }
    }

    ## Create annotated sample data after both Alamut and QCEval data loads are complete
    if ((WaitForQCLoads.wait_result.total == WaitForQCLoads.wait_result.succeeded) && (WaitForAlamutLoad.wait_result.total == WaitForAlamutLoad.wait_result.succeeded)) {
        # Wait up to 24 hours for annotation data initialization to complete
        call FASTUtils.FASTWaitForADITask {
            input:
                check_interval_secs = fast_adi_wait_interval_secs,
                max_checks = fast_adi_wait_max_intervals,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
        if (FASTWaitForADITask.num_pending >= 1) {
            call Utilities.FailTask as ADINotComplete {
                input:
                    error_message = "FAST annotation data initialization did not complete within allotted time"
            }
        }
        if (FASTWaitForADITask.num_pending <= 0) {
            # Once no variants are pending, create annotated sample data
            call FASTUtils.FASTCreateAnnotatedSampleDataTask {
                input:
                    annotated_sample_data_name = batch_name + "_" + fast_annotated_sample_data_saved_filter_name,
                    sample_data_names_and_labels = qc_fast_sample_data_name,
                    region_names_and_masks = fast_annotated_sample_data_regions,
                    scripts = fast_annotated_sample_data_scripts,
                    saved_filter_name = fast_annotated_sample_data_saved_filter_name,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    docker_image = orchutils_docker_image
            }
            # Wait for either of the annotated sample data tasks to complete
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
            # With annotation complete, export and parse results
            if (FASTWaitForAnnotatedSampleDataTask.wait_result.total == FASTWaitForAnnotatedSampleDataTask.wait_result.succeeded) {
                call FASTUtils.FASTExportAnnotatedSampleDataTask {
                    input:
                        annotated_sample_data_name = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output,
                        format = "TXT",
                        output_basename = FASTCreateAnnotatedSampleDataTask.annotated_sample_data_name_output + ".fastexport",
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
            }
        }
    }
    
    output {
        # Filtered VCFs
        Array[File] prepped_sample_vcfs = PrepSampleVCFTask.output_vcf_gz
        # Joint VCFs
        File merged_vcf_gz = MergedVCF.output_vcf_gz
        File? collective_vcf_gz = CollectiveVCF.output_vcf_gz
        # PGx output
        Array[File]? pgx_CPIC_report = PGxWorkflowAlias.CPIC_report
        Array[File]? pgx_FDA_report = PGxWorkflowAlias.FDA_report
        Array[File]? pgx_genotype_xlsx = PGxWorkflowAlias.genotype_xlsx
        Array[File]? pgx_genotype_txt = PGxWorkflowAlias.genotype_txt
        # Risk alleles output
        Array[File]? risk_alleles_report = RiskAllelesWorkflowAlias.risk_report
        Array[File]? risk_alleles_genotype_xlsx = RiskAllelesWorkflowAlias.genotype_xlsx
        Array[File]? risk_alleles_genotype_txt = RiskAllelesWorkflowAlias.genotype_txt
        # Annotated VCFs
        File alamut_vcf_gz = AlamutBatchTask.output_vcf_gz
        Array[File] qceval_vcf_gz = QCEvalTask.output_vcf_gz
        # FAST export file
        File? fast_export_file = FASTExportAnnotatedSampleDataTask.output_file
        # FAST Parsed output and NVA report
        File? fast_parsed_output = FASTOutputParserTask.parsed_report
    }
}

task PrepSampleVCFTask {
    input {
        File input_vcf
        File input_vcf_idx
        File target_roi_bed
        File ref_fasta
        File ref_fasta_index
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String docker_image
        Int disk_size = 5 + ceil(size(input_vcf, "GB") * 2.5) + ceil(size(ref_fasta, "GB")) + ceil(size(ref_fasta_index, "GB")) + ceil(size(target_roi_bed, "GB"))
    }

    command <<<
        set -euxo pipefail

        # Filter to target roi and normalize
        bcftools view --no-version --output-type z --regions-file "~{target_roi_bed}" "~{input_vcf}" > "~{output_basename}.filtered.vcf.gz"
        bcftools norm --no-version --multiallelics - --fasta-ref "~{ref_fasta}" --output-type z "~{output_basename}.filtered.vcf.gz" \
            > "~{output_basename}.filtered.norm.vcf.gz"
        # Create index file for the filtered and normalized VCF
        bcftools index --tbi "~{output_basename}.filtered.norm.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf_gz = "~{output_basename}.filtered.norm.vcf.gz"
        File output_vcf_idx = "~{output_basename}.filtered.norm.vcf.gz.tbi"
    }    
}