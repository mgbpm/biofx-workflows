version 1.0

import "../main_wdls/MakeModelWorkflow.wdl"
import "RunSinglePRSWorkflow.wdl"
import "RawScoreWorkflow.wdl"
import "MixScoreWorkflow.wdl"
import "AdjustScoreWorkflow.wdl"
import "../../../steps/Utilities.wdl"

workflow RunPrsMixWorkflow {
    input {
        File         query_vcf
        File?        score_weights
        Array[File]  adjustment_model_manifests
        String       condition_code
        Boolean      mix_before_adjustment
        Boolean      norename              = false
        String       ubuntu_docker_image   = "ubuntu:latest"
    }

    # OPTION 1: Create a model for each weights file and adjust scores before mixing
    if (!mix_before_adjustment) {
        if (!defined(score_weights)) {
            call Utilities.FailTask as MixInputsFail {
                input:
                    error_message = "Please input a score weights file for PRSmix scoring."
            }
        }

        if (defined(score_weights)) {
            scatter (manifest in adjustment_model_manifests) {
                AdjustmentModelData single_model_data = read_json(manifest)
                String weight_filename = basename(single_model_data.variant_weights[0], ".weights.tsv||.cleaned.weights.tsv")

                call RunSinglePRSWorkflow.RunSinglePrsWorkflow as ScoreSingleManifests {
                    input:
                        query_vcf = query_vcf,
                        adjustment_model_manifest = manifest,
                        norename = true,
                        ubuntu_docker_image = ubuntu_docker_image
                }
            }

            call MixScoreWorkflow.MixScoreWorkflow as MixAdjustedScores {
                input:
                    input_scores = ScoreSingleManifests.adjusted_score,
                    score_weights = select_first([score_weights]),
                    output_basename = condition_code,
                    ubuntu_docker_image = ubuntu_docker_image
            }
        }
    }

    # OPTION 2: Create a mix model and mix scores before adjusting
    if (mix_before_adjustment) {
        if (length(adjustment_model_manifests) > 1) {
            call Utilities.FailTask as MixManifestFail {
                input:
                    error_message = "Only input one PRSmix model to mix scores before adjusting."
            }
        }

        AdjustmentModelData mix_model_data = read_json(adjustment_model_manifests[0])

        call RawScoreWorkflow.RawScoreWorkflow as ScoreMixModel {
            input:
                input_vcf = query_vcf,
                adjustment_model_manifest = adjustment_model_manifests[0],
                norename = false
        }

        call MixScoreWorkflow.MixScoreWorkflow as MixRawScores {
            input:
                input_scores = ScoreMixModel.raw_scores,
                score_weights = select_first([mix_model_data.score_weights]),
                output_basename = condition_code,
                ubuntu_docker_image = ubuntu_docker_image
        }
        
        call AdjustScoreWorkflow.AdjustScoreWorkflow as AdjustMixScores {
            input:
                output_basename = condition_code,
                input_vcf = ScoreMixModel.renamed_vcf,
                adjustment_model_manifest = adjustment_model_manifests[0],
                raw_scores = MixRawScores.mix_score,
                norename = true
        }
    }  

    output {
        File final_score = select_first([MixAdjustedScores.mix_score, AdjustMixScores.adjusted_scores])
    }
}
