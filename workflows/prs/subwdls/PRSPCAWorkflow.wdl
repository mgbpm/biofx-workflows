version 1.0

import "../palantir/ScoringTasks.wdl" as ScoringTasks
import "../palantir/PCATasks.wdl" as PCATasks
import "../tasks/HelperTasks.wdl" as HelperTasks
import "../../../steps/Utilities.wdl"
import "../tasks/PRSStructs.wdl"

workflow PRSPCAWorkflow {
    input {
        String output_basename
        File input_vcf
        File adjustment_model_manifest
        File? prs_raw_scores
        Boolean norename = false
        File renaming_lookup = "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv"
    }

    AdjustmentModelData model_data = read_json(adjustment_model_manifest)

    if (! norename) {
        # Clean up the query VCF
        call HelperTasks.RenameChromosomesInVcf as RenameVcf {
            input:
                vcf = input_vcf,
                rename = renaming_lookup
        }
    }

    File input_vcf_ = select_first([RenameVcf.renamed, input_vcf])

    call PCATasks.ArrayVcfToPlinkDataset as QueryBed {
        input:
            vcf = input_vcf_,
            pruning_sites = model_data.pca_variants,
            basename = output_basename,
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
            basename = output_basename + "_pca",
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
