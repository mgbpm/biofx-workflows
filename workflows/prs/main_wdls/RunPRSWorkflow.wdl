version 1.0

import "MakeModelWorkflow.wdl"
import "../../../steps/Utilities.wdl"
import "../subwdls/RawScoreWorkflow.wdl"
import "../subwdls/MixScoreWorkflow.wdl"
import "../subwdls/AdjustScoreWorkflow.wdl"

workflow RunPrsWorkflow {
    input {
        File         query_vcf
        File?        adjustment_model_manifest
        Array[File]? variant_weights
        File?        score_weights
        File?        pca_variants
        File?        reference_vcf
        String       condition_code
        Boolean      norename = false
        String       ubuntu_docker_image = "ubuntu:latest"
    }

    if (!defined(adjustment_model_manifest)) {
        if (!defined(variant_weights) || !defined(pca_variants) || !defined(reference_vcf)) {
            call Utilities.FailTask as NoInputsFail {
                input:
                    error_message = "Required files to create a model were not provided. Please check inputs to create a model or provide one."
            }
        }

        call MakeModelWorkflow.MakeModelWorkflow as MakeModel {
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
    File model = select_first([adjustment_model_manifest, MakeModel.model_manifest])
    AdjustmentModelData model_data = read_json(model)

    call RawScoreWorkflow.RawScoreWorkflow as RawScores {
        input:
            input_vcf = select_first([model_data.query_file, query_vcf]),
            adjustment_model_manifest = model,
            norename = norename
    }

    if (defined(model_data.score_weights)) {
        call MixScoreWorkflow.MixScoreWorkflow as MixScores {
            input:
                input_scores = RawScores.raw_scores,
                score_weights = select_first([model_data.score_weights]),
                output_basename = condition_code,
                ubuntu_docker_image = ubuntu_docker_image
        }
    }

    call AdjustScoreWorkflow.AdjustScoreWorkflow as AdjustScores {
        input:
            output_basename = condition_code,
            input_vcf = RawScores.renamed_vcf,
            adjustment_model_manifest = model,
            raw_scores = select_first([MixScores.mix_score, RawScores.raw_scores]),
            norename = true
    }
    

    output {
        Array[File] raw_scores                = RawScores.raw_scores
        File?       mix_score                 = MixScores.mix_score
        File        adjusted_score            = AdjustScores.adjusted_scores
        File        adjustment_model_manifest = model
    }
}

