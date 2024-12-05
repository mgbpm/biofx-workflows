version 1.0

# Adapted from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "../tasks/ScoringTasks.wdl"
import "PRSMixScoreWorkflow.wdl"

workflow PRSTrainMixModelWorkflow {
    input {
        # Scoring inputs
        String condition_name
        Array[File] var_weights
        File scoring_sites
        File reference_vcf
        File score_weights
        Int scoring_mem
        # Training model inputs
        File population_pcs
        # Docker images
        String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        String ubuntu_docker_image = "ubuntu:21.10"
        String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    }

    # Score reference VCF
    scatter (weights_file in var_weights) {
        call ScoringTasks.ScoreVcf as ScorePopulationVCF {
            input:
                vcf = reference_vcf,
                chromosome_encoding = "MT",
                sites = scoring_sites,
                weights = weights_file,
                basename = sub(basename(weights_file), ".var_weights.tsv", ""),
                base_mem = scoring_mem,
                docker_image = plink_docker_image
        }
    }

    # Calculate PRS mix score for population VCF
    call PRSMixScoreWorkflow.PRSMixScoreWorkflow as GetMixScore {
        input:
            condition_name = condition_name,
            raw_scores = ScorePopulationVCF.score,
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
        Array[File] raw_population_scores = ScorePopulationVCF.score
        File mixed_population_scores = GetMixScore.prs_mix_raw_score
    }
}
