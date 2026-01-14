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

    File query_vcf_ = select_first([RenameChromosomesInQueryVcf.renamed, query_vcf])
    Array[File] variant_weights_ = select_first([RenameChromosomesInWeights.renamed, variant_weights])
    File pca_variants_ = select_first([RenameChromosomesInPcaVariants.renamed, pca_variants])
    File reference_vcf_ = select_first([RenameChromosomesInReferenceVcf.renamed, reference_vcf])
    
    call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
        input:
            vcf = query_vcf_
    }
    call HelperTasks.GetBaseMemory as GetMemoryForReference {
        input:
            vcf = reference_vcf_
    }
    
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
    
    call HelperTasks.TrimPcaVariants as TrimVariants {
        input:
            pca_variants = pca_variants_,
            reference = ExtractReferenceVariants.ids,
            query = ExtractQueryVariants.ids,
            docker_image = ubuntu_docker_image
    }
    File kept_pca_variants = select_first([TrimVariants.kept_pca_variants, pca_variants_])
    
    call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
        input:
            vcf = reference_vcf_,
            pruning_sites = kept_pca_variants,
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
    }

    if (defined(score_weights)) {
        call MixScoreWorkflow.MixScoreWorkflow as GetMixScore {
            input:
                output_basename = condition_code,
                input_scores = ScoreReferenceVcf.score,
                score_weights = select_first([score_weights])
        }
    }

    call ScoringTasks.TrainAncestryModel as TrainModel {
        input:
            population_pcs = ReferencePCA.pcs,
            population_scores = select_first([GetMixScore.mix_score, ScoreReferenceVcf.score[0]]),
            output_basename = condition_code
    }

    call CalculatePopStats {
        input:
            score_file = TrainModel.adjusted_population_scores,
            docker_image = ubuntu_docker_image
    }

    call BundleAdjustmentModel {
        input:
            model_data = object {
                condition_code        : condition_code,
                parameters            : "" + TrainModel.fitted_params,
                scoring_inputs        : scoring_inputs_,
                principal_components  : "" + ReferencePCA.pcs,
                loadings              : "" + ReferencePCA.pc_loadings,
                meansd                : "" + ReferencePCA.mean_sd,
                variant_weights       : variant_weights_,
                score_weights         : "" + score_weights,
                pca_variants          : "" + kept_pca_variants,
                original_pca_variants : "" + pca_variants,
                query_file            : "" + query_vcf_,
                base_memory           : GetMemoryForReference.gigabytes,
                population_mean       : CalculatePopStats.mean,
                population_sd         : CalculatePopStats.sd
            },
            docker_image = ubuntu_docker_image
    }

    output {
        Array[File] raw_reference_scores      = ScoreReferenceVcf.score
        File?       mixed_reference_scores    = GetMixScore.mix_score
        File        adjusted_reference_scores = TrainModel.adjusted_population_scores
        File        model_manifest            = BundleAdjustmentModel.manifest
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

task CalculatePopStats {
    input {
        File   score_file
        String docker_image
        Int    disk_size = ceil(size(score_file, "GB")) + 10
        Int    mem_size = 2
        Int    preemptible = 1
    }
    
    command <<<
        awk '{ total += $26; count++ } END { print total/count }' "~{score_file}" > mean.txt

        awk '{ total += $26; count++; array[NR] = $26 }
         END {
                 mean = total / count;
                 for(i=1; i<=count; i++){
                         sumsq += (array[i] - mean)^2;
                 }
                 variance = sumsq / count;
                 print sqrt(variance);
         }' "~{score_file}" > sd.txt
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk ~{disk_size} HDD"
        memory: "~{mem_size} GB"
        preemptible: preemptible
    }

    output {
        Float mean = read_float("mean.txt")
        Float sd   = read_float("sd.txt")
    }
}