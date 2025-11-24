version 1.0

import "../main_wdls/MakeModelWorkflow.wdl"
import "RawScoreWorkflow.wdl"
import "AdjustScoreWorkflow.wdl"
import "../../../steps/Utilities.wdl"

workflow RunSinglePrsWorkflow {
    input {
        File         query_vcf
        File?        variant_weights
        File?        pca_variants
        File?        reference_vcf
        File?        adjustment_model_manifest
        String       condition_code
        Boolean      norename              = false
        String       ubuntu_docker_image   = "ubuntu:latest"
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
                variant_weights = select_all([select_first([variant_weights])]),
                pca_variants = select_first([pca_variants]),
                reference_vcf = select_first([reference_vcf]),
                query_vcf = query_vcf,
                norename = norename,
                ubuntu_docker_image = ubuntu_docker_image
        }
    }

    File model = select_first([MakeModel.adjustment_model_manifest, adjustment_model_manifest])
    AdjustmentModelData model_data = read_json(model)

    call RawScoreWorkflow.RawScoreWorkflow as GetScore {
        input:
            input_vcf = select_first([model_data.query_file, query_vcf]),
            adjustment_model_manifest = model,
            norename = false
    }

    call AdjustScoreWorkflow.AdjustScoreWorkflow as AdjustScore {
        input:
            output_basename = condition_code,
            input_vcf = GetScore.renamed_vcf,
            adjustment_model_manifest = model,
            raw_scores = GetScore.raw_scores[0],
            norename = true
    }

    output {
        File raw_scores      = GetScore.raw_scores[0]
        File adjusted_scores = AdjustScore.adjusted_scores
        File output_model    = model
    }
}

