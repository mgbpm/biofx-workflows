version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"
import "../../steps/Utilities.wdl"

workflow CTBHaplotypeCallerGvcf {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230828"
        # bcftools docker image
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        # subject, sample id and data location
        String subject_id
        String sample_id
        String sample_data_location
        Boolean fetch_cram = true
        Array[String] fetch_cram_filter_keys = [subject_id, sample_id]
        Array[FileMatcher]? fetch_cram_file_matchers
        Boolean fetch_bam = true
        Array[String] fetch_bam_filter_keys = [subject_id, sample_id]
        Array[FileMatcher]? fetch_bam_file_matchers
        Boolean fetch_files_verbose = false
        Int fetch_disk_size = 75
        # reference genome files
        String reference_build = "GRCh38"
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        # variant calling inputs
        File? scattered_calling_intervals_list
    }

    # Prefer a BAM file to avoid conversion
    if (fetch_bam) {
        call FileUtils.FetchFilesTask as FetchBam {
            input:
                data_location = sample_data_location,
                file_types = if defined(fetch_bam_file_matchers) then [] else ["bam"],
                recursive = false,
                file_match_keys = if defined(fetch_bam_file_matchers) then [] else fetch_bam_filter_keys,
                file_matchers = fetch_bam_file_matchers,
                verbose = fetch_files_verbose,
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
                file_types = if defined(fetch_cram_file_matchers) then [] else ["cram"],
                recursive = false,
                file_match_keys = if defined(fetch_cram_file_matchers) then [] else fetch_cram_filter_keys,
                file_matchers = fetch_cram_file_matchers,
                verbose = fetch_files_verbose,
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
    call HaplotypeCallerGvcfGATK4.HaplotypeCallerGvcf_GATK4 {
        input:
            input_bam = sample_bam,
            input_bam_index = sample_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scattered_calling_intervals_list = select_first([scattered_calling_intervals_list]),
            make_gvcf = true
    }

    output {
        # haplotype caller output
        File? vcf = HaplotypeCallerGvcf_GATK4.output_vcf
    }
}

