version 1.0

# Adapted from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "../tasks/ScoringTasks.wdl"
import "../tasks/PRSStructs.wdl"
import "../../../steps/Utilities.wdl"
import "PRSRawScoreWorkflow.wdl"
import "PRSMixScoreWorkflow.wdl"

workflow PRSTrainMixModelWorkflow {
    input {
        # Scoring inputs
        String condition_name
        Array[File] var_weights
        File scoring_sites
        File population_vcf
        File score_weights
        # Training model inputs
        File population_pcs
        # Docker images
        String python_docker_image = "python:3.9.10"
        String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        String ubuntu_docker_image = "ubuntu:21.10"
        String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    }

    AdjustmentModelData pre_model_data = object {
        parameters: "",
        training_variants: scoring_sites,
        principal_components: "",
        loadings: "",
        meansd: "",
        var_weights: var_weights,
        score_weights: score_weights
        pca_variants: "",
        original_pca_variants: "",
        base_memory: 10
    }

    # Score population VCF with each var weights file
    call PRSRawScoreWorkflow.PRSRawScoreWorkflow as ScorePopulationVCF {
        input:
            input_vcf = population_vcf,
            adjustment_model_manifest = pre_model_data,
            python_docker_image = python_docker_image,
            plink_docker_image = plink_docker_image
    }

    # Calculate PRS mix score for population VCF
    call PRSMixScoreWorkflow.PRSMixScoreWorkflow as GetMixScore {
        input:
            condition_name = condition_name,
            raw_scores = ScorePopulationVCF.prs_raw_scores,
            score_weights = score_weights,
            ubuntu_docker_image = ubuntu_docker_image
    }

    # Train the ancestry adjustment model
    call ScoringTasks.TrainAncestryModel {
        input:
            population_pcs = population_pcs,
            population_scores = GetMixScore.prs_mix_raw_score,
            output_basename = condition_name + "_MixModel",
            docker_image = tidyverse_docker_image
    }

    output {
        # Model outputs
        File fitted_params = TrainAncestryModel.fitted_params
        Array[File] sites_used_in_scoring = ScorePopulationVCF.sites_scored
        File adjusted_population_scores = TrainAncestryModel.adjusted_population_scores
        Boolean fit_converged = TrainAncestryModel.fit_converged
        # Scoring outputs
        Array[File] raw_population_scores = ScorePopulationVCF.prs_raw_scores
        File mixed_population_scores = GetMixScore.prs_mix_raw_score
    }
}
