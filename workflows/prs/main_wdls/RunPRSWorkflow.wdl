version 1.0

import "MakeModelWorkflow.wdl"
import "../subwdls/RawScoreWorkflow.wdl"
import "../subwdls/MixScoreWorkflow.wdl"
import "../subwdls/AdjustScoreWorkflow.wdl"
import "../../../steps/Utilities.wdl"

workflow RunPrsWorkflow {
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

    # OPTION 1: Score a single weights model
    if (length(select_first([variant_weights])) == 1) {
        if (!defined(adjustment_model_manifests)) {
            call MakeModelWorkflow.MakeModelWorkflow as MakeSingleWeightsModel {
                input:
                    condition_code = condition_code,
                    variant_weights = select_first([variant_weights]),
                    pca_variants = select_first([pca_variants]),
                    reference_vcf = select_first([reference_vcf]),
                    query_vcf = query_vcf,
                    norename = norename,
                    ubuntu_docker_image = ubuntu_docker_image
            }
        }

        File? single_weights_model = select_first([
            MakeSingleWeightsModel.adjustment_model_manifest,
            adjustment_model_manifests])[0]
        AdjustmentModelData single_weights_model_data = read_json(select_first([single_weights_model]))

        call RawScoreWorkflow.RawScoreWorkflow as ScoreSingleWeightsModel {
            input:
                input_vcf = select_first([single_weights_model_data.query_file, query_vcf]),
                adjustment_model_manifest = select_first([single_weights_model]),
                norename = false
        }

        call AdjustScoreWorkflow.AdjustScoreWorkflow as AdjustSingleWeightsModel {
            input:
                output_basename = condition_code,
                input_vcf = ScoreSingleWeightsModel.renamed_vcf,
                adjustment_model_manifest = select_first([single_weights_model]),
                raw_scores = ScoreSingleWeightsModel.raw_scores[0],
                norename = true
        }
    }


    # OPTION 2: Score a mix model
    if (length(select_first([variant_weights])) > 1) {
        if (!defined(score_weights) && !mix_before_adjustment) {
            call Utilities.FailTask as MixInputsFail {
                input:
                    error_message = "Please input a score weights file for PRSmix scoring."
            }
        }

        # OPTION 2A: Create a model for each weights file and adjust scores before mixing
        if (defined(score_weights) && !mix_before_adjustment) {
            if (!defined(adjustment_model_manifests)) {
                call MakeModelWorkflow.MakeModelWorkflow as MakeAllWeightsModels {
                    input:
                        condition_code = condition_code,
                        variant_weights = select_first([variant_weights]),
                        pca_variants = select_first([pca_variants]),
                        reference_vcf = select_first([reference_vcf]),
                        query_vcf = query_vcf,
                        norename = norename,
                        ubuntu_docker_image = ubuntu_docker_image
                }
            }

            Array[File?] all_weights_models = select_first([
                MakeAllWeightsModels.adjustment_model_manifest,
                select_first([adjustment_model_manifests])
            ])

            scatter (raw_mix_model in all_weights_models) {
                AdjustmentModelData raw_mix_model_data = read_json(select_first([raw_mix_model]))

                call RawScoreWorkflow.RawScoreWorkflow as ScoreAllWeightsModel {
                    input:
                        input_vcf = select_first([raw_mix_model_data.query_file, query_vcf]),
                        adjustment_model_manifest = select_first([raw_mix_model]),
                        norename = false
                }

                call AdjustScoreWorkflow.AdjustScoreWorkflow as AdjustAllWeightsModel {
                    input:
                        output_basename = condition_code,
                        input_vcf = ScoreAllWeightsModel.renamed_vcf,
                        adjustment_model_manifest = select_first([raw_mix_model]),
                        raw_scores = ScoreAllWeightsModel.raw_scores[0],
                        norename = true
                }
            }

            call MixScoreWorkflow.MixScoreWorkflow as MixAdjustedScores {
                input:
                    input_scores = AdjustAllWeightsModel.adjusted_scores,
                    score_weights = select_first([score_weights]),
                    output_basename = condition_code,
                    ubuntu_docker_image = ubuntu_docker_image
            }
        }

        # OPTION 2B: Create a mix model and mix scores before adjusting
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
    }    

    output {
        Array[File] raw_scores = select_first([
            ScoreSingleWeightsModel.raw_scores,
            ScoreAllWeightsModel.raw_scores,
            ScoreMixModel.raw_scores
        ])

        File? mix_score = select_first([
            MixAdjustedScores.mix_score,
            MixRawScores.mix_score
        ])

        Array[File?] adjusted_scores = select_first([
            [AdjustSingleWeightsModel.adjusted_scores],
            AdjustAllWeightsModel.adjusted_scores,
            [AdjustMixScores.adjusted_scores]
        ])

        Array[File?] output_manifests = select_first([
            [single_weights_model],
            all_weights_models,
            [mix_model]
        ])
    }
}

