version 1.0

import "../main_wdls/MakeModelWorkflow.wdl"
import "RawScoreWorkflow.wdl"
import "AdjustScoreWorkflow.wdl"
import "../../../steps/Utilities.wdl"

workflow RunSinglePrsWorkflow {
    input {
        File         query_vcf
        File         adjustment_model_manifest
        String       condition_code
        Boolean      norename              = false
        String       ubuntu_docker_image   = "ubuntu:latest"
    }

    call RawScoreWorkflow.RawScoreWorkflow as GetScore {
        input:
            input_vcf = query_vcf,
            adjustment_model_manifest = adjustment_model_manifest,
            norename = false
    }

    call AdjustScoreWorkflow.AdjustScoreWorkflow as AdjustScore {
        input:
            output_basename = condition_code,
            input_vcf = GetScore.renamed_vcf,
            adjustment_model_manifest = adjustment_model_manifest,
            raw_scores = GetScore.raw_scores[0],
            norename = true
    }

    output {
        File raw_scores      = GetScore.raw_scores[0]
        File adjusted_scores = AdjustScore.adjusted_scores
    }
}

