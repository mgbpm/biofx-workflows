version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/ScoringTasks.wdl" as ScoringTasks
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/PCATasks.wdl" as PCATasks
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/HelperTasks.wdl" as HelperTasks
import "../../../steps/Utilities.wdl"
import "../tasks/PRSStructs.wdl"

workflow PRSPCAWorkflow {
    input {
        # PCA and adjustment inputs
        String condition_name
        File input_vcf
        File adjustment_model_manifest
        File? prs_raw_scores
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

    call PCATasks.ArrayVcfToPlinkDataset as QueryBed {
        input:
            vcf = RenameVcf.renamed,
            pruning_sites = model_data.pca_variants,
            basename = condition_name,
            mem = model_data.base_memory
    }

    # Run PCA with query VCF
    call PCATasks.ProjectArray as ProjectPCA {
        input:
            bim = QueryBed.bim,
            bed = QueryBed.bed,
            fam = QueryBed.fam,
            pc_loadings = model_data.loadings,
            pc_meansd = model_data.meansd,
            basename = condition_name + "_pca",
            mem = model_data.base_memory
    }

    # Plot PCA
    call PCATasks.MakePCAPlot as PCAPlot {
        input:
            population_pcs = model_data.principal_components,
            target_pcs = ProjectPCA.projections,
            output_basename = condition_name,
            docker_image = tidyverse_docker_image
    }

    # Adjust PRS raw scores with model and PCA
    if (defined(prs_raw_scores)) {
        call ScoringTasks.AdjustScores as GetAdjustedScores {
            input:
                fitted_model_params = model_data.parameters,
                pcs = ProjectPCA.projections,
                scores = select_first([prs_raw_scores])
        }
    }
    
    output {
        File pc_projection = ProjectPCA.projections
        File pc_plot = PCAPlot.pca_plot
        File? adjusted_scores = GetAdjustedScores.adjusted_scores
    }
}