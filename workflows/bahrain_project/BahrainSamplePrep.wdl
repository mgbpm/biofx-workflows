version 1.0

import "../../steps/Utilities.wdl"
import "../../steps/FileUtils.wdl"
import "../../steps/VCFUtils.wdl"

workflow BahrainSamplePrepWorkflow {
    input {
        # Sample data inputs
        Array[String] sample_ids
        Array[String] collaborator_sample_ids
        Array[String] data_bucket
        String batch_name
        File target_roi_bed
        String pipeline_to_run
        # Miscellaneous docker images
        String python_docker_image = "python:3.10"
        String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20240625"
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Reference genome files
        File ref_fasta
        File ref_fasta_index
    }

    # Test general inputs for either pipeline
    if ((pipeline_to_run != "monogenic") && (pipeline_to_run != "screening")) {
        String PipelineNameFail = "pipeline_to_run to run should be either 'monogenic' or 'screening.'"
    }
    if (length(sample_ids) != length(collaborator_sample_ids)) {
        String CollaboratorIDsLengthFail = "The length of collaborator_sample_ids does not match the length of sample_ids."
    }
    # Fail the pipeline if inputs are incorrect
    String input_error_message = select_first([PipelineNameFail, CollaboratorIDsLengthFail, ""])
    if (input_error_message != "") {
        call Utilities.FailTask as InputParameterError {
            input:
                error_message = input_error_message
        }
    }

    # Fetch CRAM, VCF, and index files
    if (input_error_message == "") {
        scatter (i in range(length(sample_ids))) {
            call FileUtils.FetchFilesTask as FetchFiles {
                input:
                    data_location = data_bucket[i],
                    recursive = true,
                    file_types = [ "vcf" ],
                    file_match_keys = [ sample_ids[i] ],
                    docker_image = orchutils_docker_image,
                    disk_size = 50,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name
            }
            if (!defined(FetchFiles.vcf)) {
                call Utilities.FailTask as VCFFileNotFound {
                    input:
                        error_message = "VCF file for sample " + sample_ids[i] + " not found in " + data_bucket[i]
                }
            }
            if (!defined(FetchFiles.vcf_index)) {
                call Utilities.FailTask as VCFIndexNotFound {
                    input:
                        error_message = "VCF index for sample " + sample_ids[i] + " not found in " + data_bucket[i]
                }
            }
        }
    }
    # Coerce object types to Array[File] for future tasks
    Array[File] fetched_vcfs = select_all(select_first([FetchFiles.vcf]))
    Array[File] fetched_vcf_idx = select_all(select_first([FetchFiles.vcf_index]))

    # Filter vcfs, merge them, and create collective vcf
    scatter (i in range(length(fetched_vcfs))) {
        String sample_output_name = collaborator_sample_ids[i] + "_" + batch_name + "_"
        call PrepSampleVCFTask {
            input:
                input_vcf = fetched_vcfs[i],
                input_vcf_idx = fetched_vcf_idx[i],
                target_roi_bed = target_roi_bed,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                output_basename = sample_output_name,
                docker_image = bcftools_docker_image
        }
    }
    call VCFUtils.MergeVCFsTask as MergedVCF {
        input:
            input_vcfs = PrepSampleVCFTask.output_vcf_gz,
            input_vcfs_idx = PrepSampleVCFTask.output_vcf_idx,
            sorted = true,
            output_basename = batch_name + ".merged",
            docker_image = bcftools_docker_image
    }
    if (pipeline_to_run == "screening") {
        call VCFUtils.MakeCollectiveVCFTask as CollectiveVCF {
            input:
                input_vcf = MergedVCF.output_vcf_gz,
                docker_image = python_docker_image
        }
    }

    output {
        # Per sample VCFs
        Array[File] sample_vcfs = PrepSampleVCFTask.output_vcf_gz
        Array[File] sample_vcf_idx = PrepSampleVCFTask.output_vcf_idx
        # Merged VCF
        File merged_vcf = MergedVCF.output_vcf_gz
        # Collective VCF
        File? collective_vcf = CollectiveVCF.output_vcf_gz
    }
}

task PrepSampleVCFTask {
    input {
        File input_vcf
        File input_vcf_idx
        File target_roi_bed
        File ref_fasta
        File ref_fasta_index
        String output_basename
        String docker_image
        Int disk_size = 5 + ceil(size(input_vcf, "GB") * 2.5) + ceil(size(ref_fasta, "GB")) + ceil(size(ref_fasta_index, "GB")) + ceil(size(target_roi_bed, "GB"))
    }

    command <<<
        set -euxo pipefail

        mkdir WORK
        mkdir OUTPUTS

        # Filter to target roi and normalize
        bcftools view --no-version --output-type z --regions-file "~{target_roi_bed}" "~{input_vcf}" > "WORK/~{output_basename}.filtered.vcf.gz"
        bcftools norm --no-version --multiallelics - --fasta-ref "~{ref_fasta}" --output-type z "WORK/~{output_basename}.filtered.vcf.gz" \
            > "WORK/~{output_basename}.filtered.norm.vcf.gz"
        # Change all instances of "Number=G" to "Number=.""
        bcftools view --header-only "WORK/~{output_basename}.filtered.norm.vcf.gz" > "WORK/old_header.txt"
        sed 's/Number=G/Number=\./' "WORK/old_header.txt" > "WORK/new_header.txt"
        bcftools reheader "WORK/~{output_basename}.filtered.norm.vcf.gz" \
            --header "WORK/new_header.txt" > "OUTPUTS/~{output_basename}.vcf.gz"
        # Create index file for final output vcf
        bcftools index --tbi "OUTPUTS/~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf_gz = "OUTPUTS/~{output_basename}.vcf.gz"
        File output_vcf_idx = "OUTPUTS/~{output_basename}.vcf.gz.tbi"
    }    
}