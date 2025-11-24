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
        Array[File]? variant_weights
        File?        score_weights
        File?        pca_variants
        File?        reference_vcf
        Array[File]? adjustment_model_manifests
        String       condition_code
        Boolean      mix_before_adjustment = true
        Boolean      norename              = false
        String       ubuntu_docker_image   = "ubuntu:latest"
    }

    # If a model manifest is not provided, check for the other inputs
    if (!defined(adjustment_model_manifests)) {
        if (!defined(variant_weights) || !defined(pca_variants) || !defined(reference_vcf)) {
            call Utilities.FailTask as NoInputsFail {
                input:
                    error_message = "Required files to create a model were not provided. Please check inputs to create a model or provide one."
            }
        }
    }

    # OPTION 1: Create a model for each weights file and adjust scores before mixing
    if (!mix_before_adjustment) {
        if (!defined(score_weights)) {
            call Utilities.FailTask as MixInputsFail {
                input:
                    error_message = "Please input a score weights file for PRSmix scoring."
            }
        }

        if (defined(score_weights) && !mix_before_adjustment) {
            if (defined(adjustment_model_manifests)) {
                Array[File] manifests = select_first([adjustment_model_manifests])

                scatter (manifest in manifests) {
                    call RunSinglePRSWorkflow.RunSinglePrsWorkflow as ScoreSingleManifests {
                        input:
                            query_vcf = query_vcf,
                            adjustment_model_manifest = manifest,
                            condition_code = condition_code,
                            norename = true,
                            ubuntu_docker_image = ubuntu_docker_image
                    }
                }
            }

            if (!defined(adjustment_model_manifests)) {
                Array[File] weights_files = select_first([variant_weights])

                scatter (weights_file in weights_files) {
                    call RunSinglePRSWorkflow.RunSinglePrsWorkflow as ScoreSingleWeights {
                        input:
                            query_vcf = query_vcf,
                            variant_weights = weights_file,
                            pca_variants = pca_variants,
                            reference_vcf = reference_vcf,
                            condition_code = condition_code,
                            ubuntu_docker_image = ubuntu_docker_image
                    }
                }
            }
        }

        Array[File] single_adjusted_scores = select_first([
            ScoreSingleManifests.adjusted_scores,
            ScoreSingleWeights.adjusted_scores
        ])
        Array[File] single_raw_scores = select_first([
            ScoreSingleManifests.raw_scores,
            ScoreSingleWeights.raw_scores
        ])
        Array[File] single_models = select_first([
            ScoreSingleManifests.output_model,
            ScoreSingleWeights.output_model
        ])

        call MixScoreWorkflow.MixScoreWorkflow as MixAdjustedScores {
            input:
                input_scores = single_adjusted_scores,
                score_weights = select_first([score_weights]),
                output_basename = condition_code,
                ubuntu_docker_image = ubuntu_docker_image
        }
    }

    # OPTION 2: Create a mix model and mix scores before adjusting
    if (mix_before_adjustment) {
        if (!defined(adjustment_model_manifests)) {
            call MakeModelWorkflow.MakeModelWorkflow as MakeMixModel {
                input:
                    condition_code = condition_code,
                    variant_weights = select_first([variant_weights]),
                    pca_variants = select_first([pca_variants]),
                    reference_vcf = select_first([reference_vcf]),
                    query_vcf = query_vcf,
                    score_weights = select_first([score_weights]),
                    norename = norename,
                    ubuntu_docker_image = ubuntu_docker_image
            }
        }

        File? mix_model = select_first([
            MakeMixModel.adjustment_model_manifest,
            adjustment_model_manifests
        ])[0]
        AdjustmentModelData mix_model_data = read_json(select_first([mix_model]))

        call RawScoreWorkflow.RawScoreWorkflow as ScoreMixModel {
            input:
                input_vcf = select_first([mix_model_data.query_file, query_vcf]),
                adjustment_model_manifest = select_first([mix_model]),
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
                adjustment_model_manifest = select_first([mix_model]),
                raw_scores = MixRawScores.mix_score,
                norename = true
        }
    }  

    output {
        File final_score = select_first([
            MixAdjustedScores.mix_score,
            AdjustMixScores.adjusted_scores
        ])

        Array[File] final_models = select_first([
            single_models,
            select_all([mix_model]),
            adjustment_model_manifests
        ])
    }
}

