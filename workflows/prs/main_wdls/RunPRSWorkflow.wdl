version 1.0

import "MakeModelWorkflow.wdl"
import "PrepareInputsWorkflow.wdl"
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
        String?      ref_source
        String?      ref_target
        String       condition_code
        Boolean      norename            = false
        String       prs_docker_image    = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/prs:20250515"
        String       ubuntu_docker_image = "ubuntu:latest"
        String       workspace
    }

    if (!defined(adjustment_model_manifest)) {
        call PrepareInputsWorkflow.PreparePrsInputs as PrepareInputs {
            input:
                variant_weights = select_first([variant_weights]),
                pca_variants = select_first([pca_variants]),
                workspace = workspace,
                source = select_first([ref_source]),
                target = select_first([ref_target]),
                norename = norename,
                query_vcfs = [query_vcf],
                prs_docker_image = prs_docker_image
        }

        Boolean rename_done = true

        if (defined(score_weights)) {
            call MakeModelWorkflow.MakeModelWorkflow as MakeMixModel {
                input:
                    condition_code = condition_code,
                    variant_weights = PrepareInputs.renamed_variant_weights,
                    pca_variants = select_first([PrepareInputs.kept_pca_variants]),
                    reference_vcf = PrepareInputs.reference_vcf,
                    query_vcf = PrepareInputs.renamed_query_vcfs[0],
                    score_weights = select_first([score_weights]),
                    norename = rename_done,
                    ubuntu_docker_image = ubuntu_docker_image
            }
        }
        if (!defined(score_weights)) {
            call MakeModelWorkflow.MakeModelWorkflow as MakeModel {
                input:
                    condition_code = condition_code,
                    variant_weights = PrepareInputs.renamed_variant_weights,
                    pca_variants = select_first([PrepareInputs.kept_pca_variants]),
                    reference_vcf = PrepareInputs.reference_vcf,
                    query_vcf = PrepareInputs.renamed_query_vcfs[0],
                    norename = rename_done,
                    ubuntu_docker_image = ubuntu_docker_image
            }
        }
    }
    File model = select_first([adjustment_model_manifest, MakeMixModel.model_manifest, MakeModel.model_manifest])
    AdjustmentModelData model_data = read_json(model)
    Boolean norename_ = select_first([rename_done, norename])

    call RawScoreWorkflow.RawScoreWorkflow as RawScores {
        input:
            input_vcf = select_first([model_data.query_file, query_vcf]),
            adjustment_model_manifest = model,
            norename = norename_
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
            input_scores = select_first([MixScores.mix_score, RawScores.raw_scores]),
            norename = true
    }
    

    output {
        Array[File]  raw_scores                = RawScores.raw_scores
        File?        mix_score                 = MixScores.mix_score
        File         adjusted_score            = AdjustScores.adjusted_scores
        File?        model_manifest            = select_first([MakeMixModel.model_manifest, MakeModel.model_manifest, adjustment_model_manifest])
        File?        kept_pca_variants         = PrepareInputs.kept_pca_variants
        Array[File]? renamed_variant_weights   = PrepareInputs.renamed_variant_weights
        File?        reference_vcf             = PrepareInputs.reference_vcf
        Array[File]? renamed_query_vcf         = PrepareInputs.renamed_query_vcfs
    }
}

