version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/subwdls/RawScoreWorkflow.wdl" as RawScoreWorkflow
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/subwdls/MixScoreWorkflow.wdl" as MixScoreWorkflow
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/subwdls/AdjustScoreWorkflow.wdl" as AdjustScoreWorkflow
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/tasks/PRSStructs.wdl" as PRSStructs
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/tasks/HelperTasks.wdl" as HelperTasks

workflow PrsScoringWorkflow {
    input {
        File input_vcf
        Array[File] model_manifests
        Boolean norename = false
        Boolean perform_adjustment = true
        File renaming_lookup = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
        String ubuntu_docker_image = "ubuntu:latest"
    }

    scatter (i in range(length(model_manifests))) {
        AdjustmentModelData model_data = read_json(model_manifests[i])

        String condition_code = model_data.condition_code

        call RawScoreWorkflow.RawScoreWorkflow as RawScores {
            input:
                input_vcf = input_vcf,
                adjustment_model_manifest = model_manifests[i],
                norename = norename,
                renaming_lookup = renaming_lookup
        }

        if (defined(model_data.score_weights)) {
            call MixScoreWorkflow.MixScoreWorkflow as MixScores {
                input:
                    output_basename = condition_code,
                    raw_scores = RawScores.prs_raw_scores,
                    score_weights = select_first([model_data.score_weights])
            }
        }

        if (perform_adjustment) {
            call AdjustScoreWorkflow.AdjustScoreWorkflow as PerformPCA {
                input:
                    output_basename = condition_code,
                    input_vcf = input_vcf,
                    adjustment_model_manifest = model_manifests[i],
                    prs_raw_scores = select_first([MixScores.prs_mix_raw_score, RawScores.prs_raw_scores]),
                    norename = norename,
                    renaming_lookup = renaming_lookup
            }
        }
    }

    output {
        Array[Array[File]] prs_raw_scores = RawScores.prs_raw_scores
        Array[File?] prs_mix_raw_score = MixScores.prs_mix_raw_score
        Array[File?] prs_adjusted_score = PerformPCA.adjusted_scores
    }
}

