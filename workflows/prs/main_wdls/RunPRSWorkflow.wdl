version 1.0

import "MakeModelWorkflow.wdl"
import "../subwdls/RawScoreWorkflow.wdl"
import "../subwdls/MixScoreWorkflow.wdl"
import "../subwdls/AdjustScoreWorkflow.wdl"

workflow RunPrsWorkflow {
    input {
        File        query_vcf
        Array[File] variant_weights
        File?       score_weights
        File        pca_variants
        File        reference_vcf
        File?       adjustment_model_manifest
        String      condition_code
        Boolean     norename = false
        String      ubuntu_docker_image = "ubuntu:latest"
    }

    if (!defined(adjustment_model_manifest)) {
        call MakeModelWorkflow.MakeModelWorkflow as MakeModel {
            input:
                condition_code = condition_code,
                variant_weights = variant_weights,
                pca_variants = pca_variants,
                reference_vcf = reference_vcf,
                query_file = query_vcf,
                score_weights = select_first([score_weights]),
                norename = norename,
                ubuntu_docker_image = ubuntu_docker_image
        }
    }

    AdjustmentModelData model_data = read_json(select_first([MakeModel.adjustment_model_manifest, adjustment_model_manifest]))

    call RawScoreWorkflow.RawScoreWorkflow as RawScores {
        input:
            input_vcf = model_data.query_file,
            adjustment_model_manifest = MakeModel.adjustment_model_manifest,
            norename = false
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
            adjustment_model_manifest = MakeModel.adjustment_model_manifest,
            raw_scores = select_first([MixScores.mix_score, RawScores.raw_scores]),
            norename = true
    }
    

    output {
        Array[File] raw_scores = RawScores.raw_scores
        File? mix_score = MixScores.mix_score
        File? adjusted_score = AdjustScores.adjusted_scores
        File? adjustment_model_manifest = MakeModel.adjustment_model_manifest
    }
}

