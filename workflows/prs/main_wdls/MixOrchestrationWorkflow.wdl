version 1.0

import "../../../steps/Utilities.wdl"
import "../../../steps/FileUtils.wdl"
import "../../lowpassimputation/Glimpse2Imputation.wdl"
import "RunPRSWorkflow.wdl"
import "../tasks/Structs.wdl"
import "../tasks/HelperTasks.wdl"

workflow MixOrchestrationWorkflow {
    input {
        String      data_location
        String      sample_id
        String      subject_id
        File        glimpse_reference_chunks
        File        ref_fasta
        File        ref_fai
        File        ref_dict
        File        gnomadAF_ref_vcf
        String      prs_test_code
        Array[File] model_manifests
        File        conditions_config
        File        percentiles
        File        renaming_lookup              = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
        String      orchutils_docker_image       = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20250203"
        String      ubuntu_docker_image          = "ubuntu:latest"
        String      glimpse_docker_image         = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String      glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        String      prs_docker_image             = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
        String      gcp_project_id               = "mgb-lmm-gcp-infrast-1651079146"
        String      workspace_name
    }

    call FileUtils.FetchFilesTask as FetchFiles {
        input:
            data_location = data_location,
            file_types = ["cram", "crai"],
            recursive = false,
            file_match_keys = [subject_id, sample_id],
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            disk_size = 75
    }

    call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
        input:
            reference_chunks = glimpse_reference_chunks,
            crams = select_all([FetchFiles.cram]),
            cram_indices = select_all([FetchFiles.crai]),
            sample_ids = [sample_id],
            fasta = ref_fasta,
            fasta_index = ref_fai,
            ref_dict = ref_dict,
            output_basename = subject_id + "_" + sample_id + "_" + prs_test_code,
            af_cutoff = ">=0.0001",
            gnomadAF_ref_vcf = gnomadAF_ref_vcf,
            impute_reference_only_variants = false,
            call_indels = false,
            keep_monomorphic_ref_sites = false,
            collect_qc_metrics = true,
            preemptible = 1,
            docker = glimpse_docker_image,
            docker_extract_num_sites_from_reference_chunk = glimpse_extract_docker_image
    }

    scatter (i in range(length(model_manifests))) {
        AdjustmentModelData model_data = read_json(model_manifests[i])

        String condition_code = model_data.condition_code

        call RunPRSWorkflow.RunPrsWorkflow as GetScores {
            input:
                query_vcf = RunGlimpse.imputed_afFiltered_vcf,
                adjustment_model_manifest = model_manifests[i],
                condition_code = condition_code,
                norename = false,
                ubuntu_docker_image = ubuntu_docker_image
        }
    }

    call SummarizeScores {
        input:
            condition_codes = condition_code,
            scores = GetScores.adjusted_score,
            conditions_config = conditions_config,
            percentiles = percentiles,
            basename = subject_id + "_" + sample_id + "_" + prs_test_code + "_results",
            docker_image = prs_docker_image
    }

    output {
        File               glimpse_vcf        = RunGlimpse.imputed_afFiltered_vcf
        File               glimpse_vcf_index  = RunGlimpse.imputed_afFiltered_vcf_index
        File?              glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[Array[File]] raw_scores         = GetScores.raw_scores
        Array[File?]       mix_scores         = GetScores.mix_score
        Array[File]        adjusted_scores    = GetScores.adjusted_score
        File               risk_summary       = SummarizeScores.risk_summary
    }
}

task SummarizeScores {
    input {
        Array[String] condition_codes
        Array[File]   scores
        File          conditions_config
        File          percentiles
        String        basename
        String        docker_image
        Int           addldisk    = 10
        Int           mem_size    = 2
        Int           preemptible = 1
    }

    Int file_size       = ceil(size(scores, "GB") + size(conditions_config, "GB"))
    Int final_disk_size = file_size + addldisk

    command <<<
        mkdir -p WORK
        mkdir -p OUTPUT

        for s in '~{sep="' '" scores}'; do
            printf -- "${s}\n" >> WORK/scores_files_list.txt
        done

        for c in '~{sep="' '" condition_codes}'; do
            printf -- "${c}\n" >> WORK/condition_codes_list.txt
        done

        summarize_mix_scores.py \
            --condition-config "~{conditions_config}" \
            --percentiles "~{percentiles}" \
            --scores-list WORK/scores_files_list.txt \
            --condition-codes WORK/condition_codes_list.txt \
            --output-file "~{basename}.tsv" \
            --output-dir OUTPUT
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File risk_summary = "OUTPUT/~{basename}.tsv"
    }
}
