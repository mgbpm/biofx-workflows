version 1.0

import "../../steps/PRSTasks.wdl"
import "../../steps/Utilities.wdl"

workflow PRSAdjustmentWorkflow {
	input {
		String condition_name
		File var_weights
		# PCA inputs
		File? pca_projections
		File? input_vcf
		File? pc_loadings
		File? pc_meansd
		File? population_pcs
		File? pruning_sites_for_pca
		File? weights_chr_encoding
		# Score inputs
		File fitted_model_params
		File prs_raw_scores
		# Docker images
		String python_docker_image = "python:3.9.10"
		String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
		String interaction_docker_image = "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
	}

	if (!defined(pca_projections) && !defined(input_vcf)) {
		call Utilities.FailTask {
            input:
                error_message = "Either a PCA projections file or all inputs to perform PCA must be provided."
        }
	}

	# Determine chromosome encoding if not provided
	if (!defined(chr_encoding)) {
		call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
			input:
				weights = var_weights,
				docker_image = python_docker_image
		}
	}

	# Run PCA if projections not provided
	if (!defined(pca_projections)) {
		call PRSTasks.ArrayVcfToPlinkDataset as GetPlinkDataset {
			input:
			vcf = input_vcf,
			pruning_sites = pruning_sites_for_pca,
			chromosome_encoding = select_first([weights_chr_encoding, ChrEncoding.chromosome_encoding]),
			basename = condition_name,
			docker_image = plink_docker_image
		}
		call PRSTasks.ProjectArray as ProjectPCA {
			input:
				bim = GetPlinkDataset.bim,
				bed = GetPlinkDataset.bed,
				fam = GetPlinkDataset.fam,
				pc_loadings = pc_loadings,
				pc_meansd = pc_meansd,
				basename = condition_name + "_pca",
				docker_image = flash_pca_docker_image
		}
		call PRSTasks.MakePCAPlot as PCAPlot {
			input:
				population_pcs = population_pcs,
				target_pcs = ProjectPCA.projections,
				docker_image = tidyverse_docker_image
		}
	}

	# Adjust score with model and PCA
	call PRSTasks.AdjustScores as GetAdjustedScores {
		input:
			fitted_model_params = fitted_model_params,
			pcs = select_first([pca_projections, ProjectPCA.projections]),
			scores = prs_raw_scores,
			docker_image = tidyverse_docker_image
	}
	
	output {
		# PCA output
		File? pc_projection = ProjectPCA.projections
		File? pc_plot = PCAPlot.pca_plot
		# Adjusted scores
		File adjusted_scores = GetAdjustedScores.adjusted_scores
	}
}