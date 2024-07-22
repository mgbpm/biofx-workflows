version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"
import "../../steps/Utilities.wdl"
#import v2.2.1
import "https://raw.githubusercontent.com/broadinstitute/warp/f0e6d797fef941c2cfea260a9dd0adcb8effe961/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl"
#import v2.1.12
#import "https://raw.githubusercontent.com/broadinstitute/warp/a4aa63170b09337a3612db6ea22e01b5b332bd54/pipelines/broad/dna_seq/germline/joint_genotyping/reblocking/ReblockGVCF.wdl"

workflow CTBHaplotypeCallerGvcf {
    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
        # bcftools docker image
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        # subject, sample id and data location
        String subject_id
        String sample_id
        String sample_data_location
        String gvcf_staging_bucket
        Boolean fetch_cram = true
        #Array[String] fetch_cram_filter_keys = [subject_id, sample_id]
        Array[String] fetch_cram_filter_keys = [subject_id]
        Array[FileMatcher]? fetch_cram_file_matchers
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

    # fetch CRAM
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

    # Validate that we got workable files from the data location
    if (!defined(FetchCram.cram)) {
        call Utilities.FailTask as MissingBamOrCramFailure {
            input:
                error_message = "CRAM file not found in ~{sample_data_location}"
        }
    }
    if (defined(FetchCram.cram) && !defined(FetchCram.crai)) {
        call Utilities.FailTask as MissingCraiFailure {
            input:
                error_message = "Index file for CRAM " + basename(select_first([FetchCram.cram])) + " not found"
        }
    }

    # Convert CRAM to BAM
    if (defined(FetchCram.cram) && defined(FetchCram.crai)) {
        call HaplotypeCallerGvcfGATK4.CramToBamTask {
            input:
                input_cram = select_first([FetchCram.cram]),
                #sample_name = basename(select_first([FetchCram.cram]), ".cram"),
                sample_name = "~{subject_id}_~{sample_id}",
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
                samtools_path = "samtools"
        }
    }

    # final inputs - either CRAM or BAM
    File sample_bam = select_first([CramToBamTask.output_bam])
    File sample_bai = select_first([CramToBamTask.output_bai])

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

    #reblock gvcf 
    call ReblockGVCF.ReblockGVCF as Reblock {
        input:
            gvcf = HaplotypeCallerGvcf_GATK4.output_vcf,
            gvcf_index = HaplotypeCallerGvcf_GATK4.output_vcf_index,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            gvcf_file_extension = ".g.vcf.gz",
            cloud_provider = "gcp",
    }


    #Transfer gvcf to staging bucket
    call FileUtils.CopyFilesTask as CopyGVCFToBucket {
    input:
        source_location = Reblock.output_vcf,
        file_types = ["g.vcf.gz"],
        file_match_keys = [],
        file_matchers = [],
        target_location = gvcf_staging_bucket,
        flatten = false,
        recursive = true,
        verbose = true,
        docker_image = orchutils_docker_image,
        disk_size = 75,
        gcp_project_id = gcp_project_id,
        workspace_name = workspace_name,
    }

   #Transfer gvcf index to staging bucket
    call FileUtils.CopyFilesTask as CopyGVCFIndexToBucket {
    input:
        source_location = Reblock.output_vcf_index,
        file_types = ["g.vcf.gz.tbi"],
        file_match_keys = [],
        file_matchers = [],
        target_location = gvcf_staging_bucket,
        flatten = false,
        recursive = true,
        verbose = true,
        docker_image = orchutils_docker_image,
        disk_size = 75,
        gcp_project_id = gcp_project_id,
        workspace_name = workspace_name,
    }

    output {
        # Haplotypecaller outputs
        File vcf = HaplotypeCallerGvcf_GATK4.output_vcf
        File vcf_index = HaplotypeCallerGvcf_GATK4.output_vcf_index
        # ReblockGVCF outputs
        File rb_vcf = Reblock.output_vcf
        File rb_vcf_index = Reblock.output_vcf_index
    }
}

