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
        String?     reported_sex
        Boolean     skip_sex_check               = false
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
        String      glimpse_docker_image         = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String      glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        String      orchutils_docker_image       = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20250203"
        String      prs_docker_image             = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
        String      python_docker_image          = "python:3.14.2"
        String      samtools_docker_image        = "biocontainers/samtools:v1.9-4-deb_cv1"
        String      ubuntu_docker_image          = "ubuntu:latest"
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

    if (defined(reported_sex) && !skip_sex_check) {
        call FindXYRatioTask {
            input:
                input_cram = select_first([FetchFiles.cram]),
                input_crai = select_first([FetchFiles.crai]),
                docker_image = samtools_docker_image
        }
        call CompareSexTask {
            input:
                xy_cov_ratio = FindXYRatioTask.xy_cov_ratio,
                input_sex = select_first([reported_sex]),
                docker_image = python_docker_image
        }
    }

    if (!defined(reported_sex) || skip_sex_check || defined(CompareSexTask.sex_guess)) {
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
    }

    scatter (i in range(length(model_manifests))) {
        AdjustmentModelData model_data = read_json(model_manifests[i])

        String condition_code = model_data.condition_code

        call RunPRSWorkflow.RunPrsWorkflow as GetScores {
            input:
                query_vcf = select_first([RunGlimpse.imputed_afFiltered_vcf]),
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
        Float?             xy_cov_ratio       = FindXYRatioTask.xy_cov_ratio
        File?              sex_guess          = CompareSexTask.sex_guess
        File?              glimpse_vcf        = RunGlimpse.imputed_afFiltered_vcf
        File?              glimpse_vcf_index  = RunGlimpse.imputed_afFiltered_vcf_index
        File?              glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[Array[File]] raw_scores         = GetScores.raw_scores
        Array[File?]       mix_scores         = GetScores.mix_score
        Array[File]        adjusted_scores    = GetScores.adjusted_score
        File               risk_summary       = SummarizeScores.risk_summary
    }
}

task FindXYRatioTask {
    input {
        File   input_cram
        File   input_crai
        String docker_image
        Int    addldisk    = 10
        Int    mem_size    = 16
        Int    preemptible = 2
    }

    Int cram_size       = ceil(size(input_cram, "GB") + size(input_crai, "GB"))
    Int final_disk_size = cram_size + addldisk

    command <<<
        set -euxo pipefail

        x_mapped_reads=$(samtools idxstats ~{input_cram} | grep -v "X_" | grep "X" | cut -f3)
        x_seq_len=$(samtools idxstats ~{input_cram} | grep -v "X_" | grep "X" | cut -f2)
        xcov=$(awk -v x=${x_mapped_reads} -v y=${x_seq_len} 'BEGIN {print x/y}')

        y_mapped_reads=$(samtools idxstats ~{input_cram} | grep -v "Y_" | grep "Y" | cut -f3)
        y_seq_len=$(samtools idxstats ~{input_cram} | grep -v "Y_" | grep "Y" | cut -f2)
        ycov=$(awk -v x=${y_mapped_reads} -v y=${y_seq_len} 'BEGIN {print x/y}')

        awk -v x=${xcov} -v y=${ycov} 'BEGIN {print x/y}' > ratio.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        Float xy_cov_ratio = read_float("ratio.txt")
    }
}

task CompareSexTask {
    input {
        Float  xy_cov_ratio
        String input_sex
        String docker_image
        Int    disk_size   = 10
        Int    mem_size    = 2
        Int    preemptible = 1
    }

    command <<<

        python -c 'with open("sex_guess.txt", "w") as f:
            if (~{xy_cov_ratio} >= 1) and (~{xy_cov_ratio} <= 4):
                sex_guess = "M"
            elif (~{xy_cov_ratio} >= 10):
                sex_guess = "F"
            else:
                sex_guess = "Undetermined"

            f.write(sex_guess)
            
            if sex_guess != "~{input_sex}":
                raise ValueError("The input sex ({1}) does not match the determined sex ({2})".format("~{input_sex}", sex_guess))
        '

    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        String sex_guess = read_string("sex_guess.txt")
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
