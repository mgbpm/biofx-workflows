version 1.0

# Adapted from commit faa824e0322e2ab455ef1cb88bc82c47d753338c of https://github.com/broadinstitute/palantir-workflows

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/ScoringTasks.wdl" as ScoringTasks
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
    }

    # Score reference VCF
    scatter (weights_file in var_weights) {
        call ScoringTasks.ScoreVcf as ScorePopulationVCF {
            input:
                vcf = reference_vcf,
                basename = sub(basename(weights_file)),
                weights = weights_file,
                base_mem = scoring_mem,
                sites = scoring_sites,
                chromosome_encoding = "MT"
        }
    }

    # Calculate PRS mix score for population VCF
    call PRSMixScoreWorkflow.PRSMixScoreWorkflow as GetMixScore {
        input:
            condition_name = condition_name,
            raw_scores = ScorePopulationVCF.score,
            score_weights = score_weights,
    }

    # Train the ancestry adjustment model
    call ScoringTasks.TrainAncestryModel {
        input:
            population_pcs = population_pcs,
            population_scores = GetMixScore.prs_mix_raw_score,
            output_basename = condition_name + "_MixModel"
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
