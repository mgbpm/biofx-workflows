version 1.0

import "../tasks/ScoringTasks.wdl"
import "../tasks/PCATasks.wdl"
import "../tasks/HelperTasks.wdl"
import "../../steps/Utilities.wdl"

workflow PRSPCAWorkflow {
    input {
        # PCA inputs
        String condition_name
        File input_vcf
        File adjustment_model_manifest
        File? var_weight_file
        String? weights_chr_encoding
        # Docker images
        String python_docker_image = "python:3.9.10"
        String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
        String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
        String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    }

    AdjustmentModelData model_data = read_json(adjustment_model_manifest)

    # Clean up the query VCF
    call HelperTasks.RenameChromosomesInVcf as RenameVcf {
        input:
            vcf = input_vcf
    }

    if (defined(model_data.query_regions)) {
        call HelperTasks.SubsetVcf as SubsetQueryVcf {
            input:
                inputvcf = RenameVcf.renamed,
                regions = select_first([model_data.query_regions]),
                label = "query"
        }
    }

    File resolved_query_vcf = select_first([SubsetQueryVcf.result, RenameVcf.renamed])

    Int base_memory = select_first([GetBaseMemoryFromVcf.gigabytes, model_data.base_memory])

    if (!defined(weights_chr_encoding)) {
        if (!defined(var_weight_file)) {
            call Utilities.FailTask as VarWeightsFail {
                input:
                    error_message = "Must have variant weights file for determining chromosome encoding."
            }
        }
        if (defined(var_weight_file)) {
            call ScoringTasks.DetermineChromosomeEncoding as ChrEncoding {
                input:
                    weights = select_first([var_weight_file]),
                    docker_image = python_docker_image
            }
        }
    }

    call PCATasks.ArrayVcfToPlinkDataset as QueryBed {
        input:
            vcf = resolved_query_vcf,
            pruning_sites = model_data.pca_variants,
            chromosome_encoding = select_first([weights_chr_encoding, ChrEncoding.chromosome_encoding]),
            basename = condition_name,
            docker_image = plink_docker_image,
            mem_size = base_memory
    }
    call PCATasks.ProjectArray as ProjectPCA {
        input:
            bim = QueryBed.bim,
            bed = QueryBed.bed,
            fam = QueryBed.fam,
            pc_loadings = model_data.loadings,
            pc_meansd = modle_data.meansd,
            basename = condition_name + "_pca",
            docker_image = flash_pca_docker_image,
            mem = base_memory
    }
    call PCATasks.MakePCAPlot as PCAPlot {
        input:
            population_pcs = model_data.principal_components,
            target_pcs = ProjectPCA.projections,
            docker_image = tidyverse_docker_image
    }
    
    output {
        File pc_projection = ProjectPCA.projections
        File pc_plot = PCAPlot.pca_plot
    }
}