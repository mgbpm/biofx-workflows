version 1.0

import "../tasks/PCATasks.wdl"
import "../tasks/ScoringTasks.wdl"
import "../subwdls/PRSTrainMixModelWorkflow.wdl"
import "../tasks/HelperTasks.wdl"
    
workflow MakeMixModelWorkflow {
    input {
        String condition_name
        Array[File] var_weights
        File pca_variants
        File reference_vcf
        File query_file
        File score_weights
        # Docker images
        String python_docker_image = "python:3.9.10"
        String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        String ubuntu_docker_image = "ubuntu:21.10"
        String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    }

    # Clean up weights, pca, and reference inputs
    scatter (i in range(length(var_weights))) {
        call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInWeights {
            input:
                tsv = var_weights[i],
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

    # Extract variant IDs from reference
    call HelperTasks.GetBaseMemory as GetMemoryForReference {
        input:
            vcf = reference_vcf
    }
    call ScoringTasks.ExtractIDsPlink as ExtractReferenceVariants {
        input:
            vcf = RenameChromosomesInReferenceVcf.renamed,
            mem = GetMemoryForReference.gigabytes
    }

    # Clean up query input
    Boolean isvcf = basename(query_file) != basename(query_file, ".vcf.gz")
    if (isvcf) {
        call HelperTasks.GetBaseMemory as GetMemoryForQueryFromVcf {
            input:
                vcf = query_file
        }
        call HelperTasks.RenameChromosomesInVcf as RenameChromosomesInQueryVcf {
            input:
                vcf = query_file
        }
        call ScoringTasks.ExtractIDsPlink as ExtractQueryVariants {
            input:
                vcf = RenameChromosomesInQueryVcf.renamed,
                mem = GetMemoryForQueryFromVcf.gigabytes
        }
    }
    if (!isvcf) {
        call HelperTasks.RenameChromosomesInTsv as RenameChromosomesInQueryVariants {
            input:
                tsv = query_file,
                skipheader = false
        }
    }

    File query_variants = select_first([ExtractQueryVariants.ids, RenameChromosomesInQueryVariants.renamed])

    call HelperTasks.TrimPCAVariants {
        input:
            pca_variants = RenameChromosomesInPcaVariants.renamed,
            reference = ExtractReferenceVariants.ids,
            query = query_variants
    }

    File kept_pca_variants = select_first([TrimPCAVariants.kept_pca_variants, RenameChromosomesInPcaVariants.renamed])

    String reference_basename = basename(reference_vcf, ".vcf.gz")

    # Perform PCA
    call PCATasks.ArrayVcfToPlinkDataset as ReferenceBed {
        input:
            vcf = select_first([SubsetReferenceVcf.result, RenameChromosomesInReferenceVcf.renamed]),
            pruning_sites = kept_pca_variants,
            basename = reference_basename,
            mem_size = GetBaseMemoryForReference.gigabytes,
            subset_to_sites = query_variants
    }
    call PCATasks.PerformPCA {
        input:
            bed = ReferenceBed.bed,
            bim = ReferenceBed.bim,
            fam = ReferenceBed.fam,
            basename = reference_basename,
            mem_size = GetMemoryForReference.gigabytes
    }

    # Train ancestry adjustment model
    call PRSTrainMixModelWorkflow.PRSTrainMixModelWorkflow as TrainModel {
        input:
            condition_name = condition_name,
            var_weights = RenameChromosomesInWeights.renamed,
            scoring_sites = query_variants,
            population_vcf = RenameChromosomesInReferenceVcf.renamed,
            score_weights = score_weights,
            population_pcs = PerformPCA.pcs,
            python_docker_image = python_docker_image,
            plink_docker_image = plink_docker_image,
            ubuntu_docker_image = ubuntu_docker_image,
            tidyverse_docker_image = tidyverse_docker_image
    }

    # Bundle model outputs in JSON
    Array[String] renamed_weights = RenameChromosomesInWeights.renamed

    call BundleAdjustmentModel {
        # Convert files to string types so a VM is not used
        input:
            model_data = object {
                parameters : "" + TrainModel.fitted_params,
                training_variants : "" + TrainModel.sites_used_in_scoring,
                principal_components : "" + PerformPCA.pcs,
                loadings : "" + PerformPCA.pc_loadings,
                meansd : "" + PerformPCA.mean_sd,
                score_weights : "" + score_weights,
                var_weights : renamed_weights,
                pca_variants : "" + kept_pca_variants,
                original_pca_variants: "" + pca_variants,
                base_memory : GetMemoryForReference.gigabytes
            },
            docker_image = ubuntu_docker_image
    }

    output {
        # Model outputs
        File adjustment_model_manifest = BundleAdjustmentModel.manifest
        Boolean fit_converged = TrainModel.fit_converged
        # Score outputs
        Array[File] raw_reference_scores = TrainModel.raw_population_scores
        File mixed_reference_scores = TrainModel.mixed_population_scores
        File adjusted_reference_scores = TrainModel.adjusted_population_scores
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
