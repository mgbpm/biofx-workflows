version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/main_wdls/MakeMixModelWorkflow.wdl" as MixModelWorkflow
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/workflows/prs/subwdls/MakeAdjustmentModelWorkflow.wdl" as SingleModelWorkflow
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/steps/Utilities.wdl" as Utilities

workflow PrsModel {
    input {
        String condition_name
        Array[File] var_weights
        File pca_variants
        File reference_vcf
        File query_file
        File? score_weights
        Boolean norename = false
    }

    # Create mix model if multiple variant weights are provided and score weights are defined
    if (length(var_weights) > 1 && defined(score_weights)) {
        call MixModelWorkflow.MakeMixModelWorkflow as MakeMixModel {
            input:
                condition_code = condition_name,
                var_weights = var_weights,
                pca_variants = pca_variants,
                reference_vcf = reference_vcf,
                query_file = query_file,
                score_weights = select_first([score_weights]),
                norename = norename
        }
    }


    # Create single score model if only one variant weight is provided and score weights are not defined
    if (length(var_weights) == 1 && !defined(score_weights)) {
        call SingleModelWorkflow.MakeAdjustmentModel as MakeModel {
            input:
                weights = var_weights,
                pca_variants = pca_variants,
                reference_vcf = reference_vcf,
                query_file = query_file,
                name = condition_name,
                norename = norename,
        }
    }

    # Fail if there is a mismatch in the number of variant weights and score weights
    if (length(var_weights) > 1 && !defined(score_weights)) {
        call Utilities.FailTask as NoScoreWeights {
            input:
                error_message = "Multiple variant weight files provided, but no score weights file specified. Please provide a score weights file for mix model scoring."
        }
    }
    if (length(var_weights) == 1 && defined(score_weights)) {
        call Utilities.FailTask as SingleWeights {
            input:
                error_message = "Single variant weight file provided, but score weights file specified. Please remove the score weights file for single model scoring or include multiple variant weights files for mix scoring."
        }
    }

    output {
        File model_manifest = select_first([MakeMixModel.adjustment_model_manifest, MakeModel.adjustment_model_manifest])
    }
}
