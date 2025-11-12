version 1.0

import "../tasks/PCATasks.wdl"
import "../tasks/ScoringTasks.wdl"
import "../subwdls/MixScoreWorkflow.wdl"
import "../tasks/HelperTasks.wdl"
import "../../../steps/Utilities.wdl"

workflow MakeModelWorkflow {
  input {
    String      condition_code
    Array[File] variant_weights
    File        pca_variants
    File        reference_vcf
    File        query_file
    File?       score_weights
    Boolean     norename            = false
    String      ubuntu_docker_image = "ubuntu:21.10"
  }

  if ((length(variant_weights) > 1) && !defined(score_weights)) {
    call Utilities.FailTask as MissingScoreWeights {
      input:
        error_message = "Missing score_weights file. It is required when multiple variant weights are specified."
    }
  }

  if ((length(variant_weights) == 1) || defined(score_weights)) {
    if (defined(score_weights)) {
      call HelperTasks.CheckInputWeightFiles {
        input:
          score_weights = select_first([score_weights]),
          variant_weights = variant_weights,
          docker_image = ubuntu_docker_image
      }
    }
  
    if (!norename) {
      scatter (i in range(length(variant_weights))) {
        call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
          input:
            tsv = variant_weights[i],
            skipheader = true
        }
      }
  
      call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInPcaVariants {
        input:
          tsv = pca_variants,
          skipheader = false
      }
  
      call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInReferenceVcf {
        input:
          vcf = reference_vcf
      }
    }
  
    Array[File] variant_weights_ = select_first([RenameChromosomesInWeights.renamed, variant_weights])
    File pca_variants_ = select_first([RenameChromosomesInPcaVariants.renamed, pca_variants])
    File reference_vcf_ = select_first([RenameChromosomesInReferenceVcf.renamed, reference_vcf])
  
    call HelperTasks.GetBaseMemory as GetMemoryForReference {
      input:
        vcf = reference_vcf
    }
  
    call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
      input:
        vcf = reference_vcf_,
        mem_size = GetMemoryForReference.gigabytes
    }
  
    Boolean isvcf = basename(query_file) != basename(query_file, ".vcf.gz")
  
    if (isvcf) {
      call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
        input:
          vcf = query_file
      }
  
      if (!norename) {
        call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
          input:
            vcf = query_file
        }
      }
  
      File query_vcf = select_first([RenameChromosomesInQueryVcf.renamed, query_file])
  
      call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
        input:
          vcf = select_first([RenameChromosomesInQueryVcf.renamed, query_file]),
          mem_size = GetMemoryForQueryFromVcf.gigabytes
      }
    }
    if (! isvcf) {
      if (!norename) {
        call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInQueryVariants {
          input:
            tsv = query_file,
            skipheader = false
        }
      }
      File query_file_ = select_first([RenameChromosomesInQueryVariants.renamed, query_file])
    }
  
    File query_variants = select_first([ExtractQueryVariants.ids, query_file_])
  
    call HelperTasks.TrimPcaVariants as TrimVariants {
      input:
        pca_variants = pca_variants_,
        reference = ExtractReferenceVariants.ids,
        query = query_variants,
        docker_image = ubuntu_docker_image
    }

  
    String reference_basename = basename(reference_vcf, ".vcf.gz")
  
    call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
      input:
        vcf = reference_vcf_,
        pruning_sites = select_first([TrimVariants.kept_pca_variants, pca_variants_]),
        subset_to_sites = query_variants,
        basename = reference_basename,
        mem_size = GetMemoryForReference.gigabytes
    }
    call PCATasks.PerformPCA {
      input:
        bed = ReferenceBed.bed,
        bim = ReferenceBed.bim,
        fam = ReferenceBed.fam,
        basename = reference_basename,
        mem_size = GetMemoryForReference.gigabytes
    }
  
    # Train the model
    scatter (weights_file in variant_weights_) {
      call ScoringTasks.ScoreVcf as ScorePopulationVcf {
        input:
          vcf = reference_vcf_,
          basename = basename(weights_file, ".tsv"),
          weights = weights_file,
          mem_size = GetMemoryForReference.gigabytes,
          sites = query_variants,
          chromosome_encoding = "MT"
      }

      Object scoring_inputs_ = object{
        variant_weights : "" + weights_file,
        training_variants: "" + ScorePopulationVcf.sites_scored
      }
    }

    if (length(variant_weights_) > 1) {
      call MixScoreWorkflow.MixScoreWorkflow as GetMixScore {
        input:
          output_basename = condition_code + "_" + basename(reference_vcf_),
          input_scores = ScorePopulationVcf.score,
          score_weights = select_first([score_weights])
      }
    }

    File population_scores = select_first([GetMixScore.mix_score, ScorePopulationVcf.score[0]])

    call ScoringTasks.TrainAncestryModel {
      input:
        population_pcs = PerformPCA.pcs,
        population_scores = population_scores,
        output_basename = condition_code
    }
  
    if (defined(score_weights)) {
      call BundleAdjustmentModel as BundleMixModel {
        input:
          model_data = object {
            condition_code : condition_code,
            parameters : "" + TrainAncestryModel.fitted_params,
            scoring_inputs : scoring_inputs_,
            principal_components : "" + PerformPCA.pcs,
            loadings : "" + PerformPCA.pc_loadings,
            meansd : "" + PerformPCA.mean_sd,
            variant_weights : "" + variant_weights_,
            score_weights : "" + score_weights,
            pca_variants : "" + TrimVariants.kept_pca_variants,
            original_pca_variants : "" + pca_variants,
            query_file : "" + select_first([query_vcf, query_file_]),
            base_memory : GetMemoryForReference.gigabytes
          },
          docker_image = ubuntu_docker_image
      }
    }
  
    if (!defined(score_weights)) {
      call BundleAdjustmentModel as BundleModel {
        input:
          model_data = object {
            condition_code : condition_code,
            parameters : "" + TrainAncestryModel.fitted_params,
            scoring_inputs : scoring_inputs_,
            principal_components : "" + PerformPCA.pcs,
            loadings : "" + PerformPCA.pc_loadings,
            meansd : "" + PerformPCA.mean_sd,
            variant_weights : "" + variant_weights_,
            pca_variants : "" + TrimVariants.kept_pca_variants,
            original_pca_variants : "" + pca_variants,
            query_file : "" + select_first([query_vcf, query_file_]),
            base_memory : GetMemoryForReference.gigabytes
          },
          docker_image = ubuntu_docker_image
      }
    }
  }

  output {
    File         adjustment_model_manifest = select_first([BundleMixModel.manifest, BundleModel.manifest])
    Array[File]? raw_reference_scores      = ScorePopulationVcf.score
    File?        mixed_reference_scores    = GetMixScore.mix_score
    File?        adjusted_reference_scores = TrainAncestryModel.adjusted_population_scores
  }
}

task BundleAdjustmentModel {
  input {
    Object model_data
    String docker_image
  }

  command {}

  output {
    File manifest = write_json(model_data)
  }

  runtime {
    docker: "~{docker_image}"
  }
}