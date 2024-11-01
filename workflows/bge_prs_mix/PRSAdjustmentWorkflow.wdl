version 1.0

import "PRSTasks.wdl"
import "../../steps/Utilities.wdl"

workflow PRSAdjustmentWorkflow {
	input {
		# Zip with condition-specific files
		File condition_file
		# PCA inputs
		File? pca_projections
		File? input_vcf
		File? pruning_sites_for_pca
		String? chr_encoding
		# PRS Mix raw scores
		File prs_mix_raw_score
		# Docker images
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

	String condition_name = sub(basename(condition_file), "\\.(tar|TAR|tar.gz|TAR.GZ)$", "")

	# Determine chromosome encoding if not provided
	if (!defined(chr_encoding)) {
		call PRSTasks.DetermineChromosomeEncoding as ChrEncoding {
			input:
				condition_zip_file = condition_file,
				docker_image = interaction_docker_image
		}
	}

	# Run PCA if projections not provided
	if (!defined(pca_projections)) {
		call PRSTasks.ArrayVCFToPlinkDataset as VCFToPlinkDataset {
			input:
			input_vcf = select_first([input_vcf]),
			pruning_sites = select_first([pruning_sites_for_pca]),
			chromosome_encoding = select_first([chr_encoding, ChrEncoding.chr_encoding]),
			output_basename = condition_name,
			docker_image = plink_docker_image
		}
		call PRSTasks.ProjectArray as ProjectPC {
			input:
				condition_zip_file = condition_file,
				bed = VCFToPlinkDataset.bed,
				bim = VCFToPlinkDataset.bim,
				fam = VCFToPlinkDataset.fam,
				output_basename = condition_name + "_pca",
				docker_image = flash_pca_docker_image
		}
		call PRSTasks.MakePCAPlot as PCAPlot {
			input:
				condition_zip_file = condition_file,
				target_pcs = ProjectPC.projections,
				output_basename = condition_name,
				docker_image = tidyverse_docker_image
		}
	}

	# Adjust score with model and PCA
	call PRSTasks.AdjustScores as GetAdjustedScores {
		input:
			condition_zip_file = condition_file,
			pcs = select_first([pca_projections, ProjectPC.projections]),
			scores = prs_mix_raw_score,
			output_basename = condition_name,
			docker_image = tidyverse_docker_image
	}
	
	output {
		# PCA output
		File? pc_projection = ProjectPC.projections
		File? pc_plot = PCAPlot.pca_plot
		# Adjusted scores
		File adjusted_scores = GetAdjustedScores.adjusted_scores
	}
}