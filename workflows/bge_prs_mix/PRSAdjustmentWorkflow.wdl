version 1.0

import "../../steps/PRSTasks.wdl"
import "../../steps/Utilities.wdl"

workflow PRSAdjustmentWorkflow {
	input {
		# Score adjustment inputs
		String condition_name
		File? var_weight_file
		File fitted_model_params
		File prs_raw_scores
		# Chromosome encoding
		String? weights_chr_encoding
		# PCA inputs
		File? pca_projections
		File? input_vcf
		File? pc_loadings
		File? pc_meansd
		File? population_pcs
		File? pruning_sites_for_pca
		# Docker images
		String python_docker_image = "python:3.9.10"
		String tidyverse_docker_image = "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
		String plink_docker_image = "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
		String flash_pca_docker_image = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
	}

	if (!defined(pca_projections) && !defined(input_vcf)) {
		call Utilities.FailTask as PCAInputFail {
            input:
                error_message = "Either a PCA projections file or all inputs to perform PCA must be provided."
        }
	}

	# Determine chromosome encoding if not provided
	if (!defined(weights_chr_encoding)) {
		if (!defined(var_weight_file)) {
			call Utilities.FailTask as VarWeightsFail {
				input:
					error_message = "Must have variant weights file for determining chromosome encoding."
			}
		}
		call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
			input:
				weights = select_first([var_weight_file]),
				docker_image = python_docker_image
		}
	}

	# Run PCA if projections not provided
	if (!defined(pca_projections)) {
		call PRSTasks.ArrayVcfToPlinkDataset as GetPlinkDataset {
			input:
			vcf = select_first([input_vcf]),
			pruning_sites = select_first([pruning_sites_for_pca]),
			chromosome_encoding = select_first([weights_chr_encoding, ChrEncoding.chromosome_encoding]),
			basename = condition_name,
			docker_image = plink_docker_image
		}
		call PRSTasks.ProjectArray as ProjectPCA {
			input:
				bim = GetPlinkDataset.bim,
				bed = GetPlinkDataset.bed,
				fam = GetPlinkDataset.fam,
				pc_loadings = select_first([pc_loadings]),
				pc_meansd = select_first([pc_meansd]),
				basename = condition_name + "_pca",
				docker_image = flash_pca_docker_image
		}
		call PRSTasks.MakePCAPlot as PCAPlot {
			input:
				population_pcs = select_first([population_pcs]),
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