version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/ScoringTasks.wdl" as ScoringTasks
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/PCATasks.wdl" as PCATasks
import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/feature/prs/gb083__TEMP__241122F195830/workflows/prs/HelperTasks.wdl" as HelperTasks
import "../../../steps/Utilities.wdl"
import "../tasks/PRSStructs.wdl"

workflow PRSPCAWorkflow {
    input {
        String condition_name
        File input_vcf
        File adjustment_model_manifest
        File? prs_raw_scores
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
    call ScoringTasks.MakePCAPlot as PCAPlot {
        input:
            population_pcs = model_data.principal_components,
            target_pcs = ProjectPCA.projections
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