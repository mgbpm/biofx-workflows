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
    File        query_vcf
    File?       score_weights
    Boolean     norename            = false
    String      ubuntu_docker_image = "ubuntu:21.10"
  }

  # If score weights are supplied, check that they match the variant weights
  if (defined(score_weights)) {
    if (length(variant_weights) == 1) {
      call Utilities.FailTask as VarWeightsFail {
        input:
          error_message = "If a score weights file is given, there must be more than one variant weights file as input."
      }
    }
    if (length(variant_weights) > 1) {
      call HelperTasks.CheckInputWeightFiles {
        input:
          score_weights = select_first([score_weights]),
          variant_weights = variant_weights,
          docker_image = ubuntu_docker_image
      }
    }
  }
  
  # Rename chromosomes in all inputs if necessary
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
    call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
      input:
        vcf = query_vcf
    }
  }

  # Resolve input files for easy reference
  File query_vcf_ = select_first([RenameChromosomesInQueryVcf.renamed, query_vcf])
  Array[File] variant_weights_ = select_first([RenameChromosomesInWeights.renamed, variant_weights])
  File pca_variants_ = select_first([RenameChromosomesInPcaVariants.renamed, pca_variants])
  File reference_vcf_ = select_first([RenameChromosomesInReferenceVcf.renamed, reference_vcf])
  
  # Get allocated memory from reference and query
  call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
    input:
      vcf = query_vcf_
  }
  call HelperTasks.GetBaseMemory as GetMemoryForReference {
    input:
      vcf = reference_vcf_
  }
  
  # Extract variants from reference and query
  call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
    input:
      vcf = query_vcf_,
      mem_size = GetMemoryForQueryFromVcf.gigabytes
  }
  call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
    input:
      vcf = reference_vcf_,
      mem_size = GetMemoryForReference.gigabytes
  }
  
  # Trim PCA variants file to variants found in query and reference
  call HelperTasks.TrimPcaVariants as TrimVariants {
    input:
      pca_variants = pca_variants_,
      reference = ExtractReferenceVariants.ids,
      query = ExtractQueryVariants.ids,
      docker_image = ubuntu_docker_image
  }
  
  # Make PLINK bed from reference VCF and perform PCA
  call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
    input:
      vcf = reference_vcf_,
      pruning_sites = select_first([TrimVariants.kept_pca_variants, pca_variants_]),
      subset_to_sites = ExtractQueryVariants.ids,
      basename = basename(reference_vcf, ".vcf.gz"),
      mem_size = GetMemoryForReference.gigabytes
  }
  call PCATasks.PerformPCA as ReferencePCA {
    input:
      bed = ReferenceBed.bed,
      bim = ReferenceBed.bim,
      fam = ReferenceBed.fam,
      basename = basename(reference_vcf, ".vcf.gz"),
      mem_size = GetMemoryForReference.gigabytes
  }
  
  scatter (weights_file in variant_weights_) {
    # Score the reference VCF per weights file
    call ScoringTasks.ScoreVcf as ScoreReferenceVcf {
      input:
        vcf = reference_vcf_,
        basename = basename(weights_file, ".tsv"),
        weights = weights_file,
        mem_size = GetMemoryForReference.gigabytes,
        sites = ExtractQueryVariants.ids,
        chromosome_encoding = "MT"
    }
    Object scoring_inputs_ = object{
      variant_weights : "" + weights_file,
      training_variants: "" + ScoreReferenceVcf.sites_scored
    }

    # If a score weights file is not supplied, either a single weights model is
    # desired or future scores will be mixed after adjusting with each weights
    # files' model
    if (!defined(score_weights)) {
      # Train a model per weights file
      call ScoringTasks.TrainAncestryModel as TrainModel {
        input:
          population_pcs = ReferencePCA.pcs,
          population_scores = ScoreReferenceVcf.score,
          output_basename = condition_code
      }
      # Bundle each weights files' model
      Array[File] weights_file_ = [weights_file]

      call BundleAdjustmentModel as BundleModel {
        input:
          model_data = object {
            condition_code : condition_code,
            parameters : "" + TrainModel.fitted_params,
            scoring_inputs : scoring_inputs_,
            principal_components : "" + ReferencePCA.pcs,
            loadings : "" + ReferencePCA.pc_loadings,
            meansd : "" + ReferencePCA.mean_sd,
            variant_weights : "" + weights_file_,
            pca_variants : "" + TrimVariants.kept_pca_variants,
            original_pca_variants : "" + pca_variants,
            query_file : "" + query_vcf_,
            base_memory : GetMemoryForReference.gigabytes
          },
          docker_image = ubuntu_docker_image
      }
    }
  }

  if (defined(score_weights)) {
    call MixScoreWorkflow.MixScoreWorkflow as GetMixScore {
      input:
        output_basename = condition_code + "_" + basename(reference_vcf_),
        input_scores = ScoreReferenceVcf.score,
        score_weights = select_first([score_weights])
    }
    # Train model on mix scores
    call ScoringTasks.TrainAncestryModel as TrainMixModel {
      input:
        population_pcs = ReferencePCA.pcs,
        population_scores = GetMixScore.mix_score,
        output_basename = condition_code
    }
    # Bundle mix model
    call BundleAdjustmentModel as BundleMixModel {
      input:
        model_data = object {
          condition_code : condition_code,
          parameters : "" + TrainMixModel.fitted_params,
          scoring_inputs : scoring_inputs_,
          principal_components : "" + ReferencePCA.pcs,
          loadings : "" + ReferencePCA.pc_loadings,
          meansd : "" + ReferencePCA.mean_sd,
          variant_weights : "" + variant_weights_,
          score_weights : "" + score_weights,
          pca_variants : "" + TrimVariants.kept_pca_variants,
          original_pca_variants : "" + pca_variants,
          query_file : "" + query_vcf_,
          base_memory : GetMemoryForReference.gigabytes
        },
        docker_image = ubuntu_docker_image
    }
  }

  output {
    Array[File]? raw_reference_scores        = ScoreReferenceVcf.score
    File?        mixed_reference_scores      = GetMixScore.mix_score

    Array[File?] adjusted_reference_scores   = select_first([
      TrainModel.adjusted_population_scores,
      [TrainMixModel.adjusted_population_scores]
    ])
    
    Array[File?] adjustment_model_manifest   = select_first([
      BundleModel.manifest,
      [BundleMixModel.manifest]
    ])
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