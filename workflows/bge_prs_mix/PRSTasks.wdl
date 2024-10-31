version 1.0

task CalculateMixScore {
	input {
		Array[File] raw_scores_files
		Int raw_scores_len = length(raw_scores_files)
		File score_weights
		String output_basename
		String docker_image
		Int disk_size = ceil(size(score_weights, "GB") * 2) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT
		mkdir -p OUTPUT

		# Extract all sample IDs from a raw score file
		score_file_array=('~{sep="' '" raw_scores_files}')
		sed '1d;' ${score_file_array[0]} | awk '{ print $1 }' > OUTPUT/sample_ids.txt

		# Set up score file headers
		printf "#IID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\n" > "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		while read line; do
			# Initialize sum of sample's raw scores
			score_sum=0

			# Add the raw score from each file to the sum of raw scores
			for c in '~{sep="' '" raw_scores_files}'; do
				pgs_id=$(basename $c .txt | cut -d "_" -f 1)
				score_weight=$(grep "${pgs_id}" "~{score_weights}" | cut -f 2)
				raw_score=$(grep "${line}" $c | cut -f 4)
				weighted_score=$(awk -v x=${score_weight} -v y=${raw_score} 'BEGIN {print x*y}')
				score_sum=$(awk -v x=${weighted_score} -v y=${score_sum} 'BEGIN {print x+y}')
			done

			# Get the weighted average for raw scores (the PRS mix score)
			weighted_avg=$(awk -v x=${score_sum} -v y="~{raw_scores_len}" 'BEGIN {print x/y}')

			# Print info for the sample
			printf "${line}\t0\t${weighted_avg}\t${score_sum}\n" >> "OUTPUT/~{output_basename}_prs_mix_score.sscore"

		done < OUTPUT/sample_ids.txt
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		File sample_ids_list = "OUTPUT/sample_ids.txt"
		File prs_mix_raw_score = "OUTPUT/~{output_basename}_prs_mix_score.sscore"
	}
}

task UnzipConditionFile {
	input {
		File zipped_condition_file
		String output_basename = sub(basename(zipped_condition_file), "\\.(tar|TAR|tar.gz|TAR.GZ)$", "")
		String docker_image
		Int disk_size = ceil(size(zipped_condition_file, "GB") * 2) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT

		tar -xf "${zipped_condition_file}" -C OUTPUT
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		Array[File] var_weights_files = glob("OUTPUT/*.var_weights.txt")
		File score_weights_file = "OUTPUT/~{output_basename}.score_weights.txt"
		File population_loadings = "OUTPUT/~{output_basename}.pc.loadings"
		File population_meansd = "OUTPUT/~{output_basename}.pc.meansd"
		File population_pcs = "OUTPUT/~{output_basename}.pc"
		File ancestry_adjustment_model = "OUTPUT/~{output_basename}.fitted_model_params.tsv"
	}	
}

task CategorizeScores {
	input {
		Array[File] adjusted_scores
		File sample_ids
		String docker_image
		Int disk_size = ceil(size(adjusted_scores, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p OUTPUT

		while read line; do
			# Write headers for individual's risk summary
			printf "Condition,Percentile,Total_Bins,Individuals_Bin\n">> OUTPUT/"${line}"_risk_summary.csv	

			for c in '~{sep="' '" adjusted_scores}'; do
				# Find the bin for the percentile
				percentile=$(grep "${line}" $c | cut -f 27)
				if [ $percentile -le 20 ]; then
					bin="Below Average"
				elif [ $percentile -le 75]; then
					bin="Average"
				elif [ $percentile -le 90]; then
					bin="Above Average"
				else
					bin="High"
				fi
				
				# Output individual's percentile and bin to csv
				printf "$(basename ${c} _adjusted_scores.tsv),${percentile},4,${bin}\n" >> OUTPUT/"${line}"_risk_summary.csv
			done
		done < "~{sample_ids}"
	>>>

	runtime {
  		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_size + "GB"
		preemptible: preemptible
	}

	output {
		Array[File] individual_risk_summaries = glob("OUTPUT/*.csv")
	}	
}

# Following PRS tasks adapted from:
# https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/ScoringTasks.wdl
# https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/PCATasks.wdl

# plink chromosome encoding rules: https://www.cog-genomics.org/plink/2.0/data#irreg_output
task DetermineChromosomeEncoding {
	input {
		File var_weights
		String output_basename
		String docker_image
	}

	command <<<
		python3 << "EOF"
			code = 'MT'
			with open("~{var_weights}") as weights_file:
				chroms = {s.split("\t")[0].split(":")[0] for i, s in enumerate(weights_file) if i > 0}
				if any('chr' in c for c in chroms):
					if 'chrM' in chroms:
						code = 'chrM'
					else:
						code = 'chrMT'
				elif 'M' in chroms:
					code = 'M'

			with open('~{output_basename}.chr_encode.txt', 'w') as write_code_file:
				write_code_file.write(f'{code}\n')
		EOF
	>>>

	runtime {
		docker : "~{docker_image}"
	}

	output {
		String chromosome_encoding = read_string("~{output_basename}.chr_encode.txt")
	}
}

task ScoreVCF {
	input {
		File input_vcf
		File var_weights
		String chromosome_encoding
		File sites
		String output_basename
		String docker_image
		Int base_mem = 16
		Int mem_size = base_mem + 2
		Int plink_mem = ceil(base_mem * 0.75 * 1000)
		Int disk_size =  (ceil(size(input_vcf, "GB")) * 3) + 20
	}
  
	command <<<
		set -euxo pipefail

		/plink2 --score "~{var_weights}" \
			header ignore-dup-ids list-variants no-mean-imputation \
			cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums \
			--set-all-var-ids "@:#:\\$1:\\$2" \
			-vcf "~{input_vcf}" \
			dosage=DS \
			--new-id-max-allele-len 1000 missing \
			--extract "~{sites}" \
			--out "~{output_basename}" \
			--memory "~{plink_mem}" \
			--output-chr "~{chromosome_encoding}"
	>>>

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
	}

	output {
		File score = "~{output_basename}.sscore"
		File log = "~{output_basename}.log"
		File sites_scored = "~{output_basename}.sscore.vars"
	}
}

task ExtractIDs {
	input {
		File input_vcf
		String chromosome_encoding
		String docker_image
		Int disk_size = 2 * ceil(size(input_vcf, "GB")) + 100
		Int mem_size = 8
		Int plink_mem = ceil(mem_size * 0.75 * 1000)
	}

	command <<<
		set -euxo pipefail

		/plink2 \
			--vcf "~{input_vcf}" \
			--set-all-var-ids "@:#:\\$1:\\$2" \
			--new-id-max-allele-len 1000 missing \
			--rm-dup exclude-all \
			--allow-extra-chr \
			--write-snplist allow-dups \
			--memory "~{plink_mem}" \
			--output-chr "~{chromosome_encoding}"
	>>>

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
	}

	output {
		File ids = "plink2.snplist"
	}
}

task ArrayVCFToPlinkDataset {
	input {
		File input_vcf
		File pruning_sites
		String chromosome_encoding
    	String output_basename
		String docker_image
		Int disk_size = (ceil(size(input_vcf, "GB")) * 3) + 20
		Int mem_size = 8
	}

	command <<<
		/plink2 \
			--vcf "~{input_vcf}" \
			--extract-intersect "~{pruning_sites}" \
			--allow-extra-chr \
			--set-all-var-ids "@:#:\\$1:\\$2" \
			--new-id-max-allele-len 1000 missing \
			--out "~{output_basename}" \
			--make-bed \
			--rm-dup force-first \
			--output-chr "~{chromosome_encoding}"
	>>>

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
	}

	output {
		File bed = "~{output_basename}.bed"
		File bim = "~{output_basename}.bim"
		File fam = "~{output_basename}.fam"
	}
}

task ProjectArray {
	input {
		File bim
		File bed
		File fam
		File pc_loadings
		File pc_meansd
		String output_basename
		String docker_image
		Int disk_size = 400
		Int mem_size = 8
		Int nthreads = 16
	}

	command <<<
		cp "~{bim}" "~{output_basename}.bim"
		cp "~{bed}" "~{output_basename}.bed"
		cp "~{fam}" "~{output_basename}.fam"

		cp "~{pc_loadings}" loadings.txt
		cp "~{pc_meansd}" meansd.txt

		# Check if .bim file, pc loadings, and pc meansd files have the same IDs
		# 1. extract IDs, removing first column of .bim file and first rows of the pc files
		awk '{print $2}' "~{output_basename}.bim" > bim_ids.txt
		awk '{print $1}' loadings.txt | tail -n +2 > pcloadings_ids.txt
		awk '{print $1}' meansd.txt | tail -n +2 > meansd_ids.txt

		diff bim_ids.txt pcloadings_ids.txt > diff1.txt
		diff bim_ids.txt meansd_ids.txt > diff2.txt
		diff pcloadings_ids.txt meansd_ids.txt > diff3.txt

		if [[ -s diff3.txt ]];then
			echo "PC loadings file and PC means file do not contain the same IDs. Check your input files and run again."
			exit 1
		fi

		# Check if diff files are not empty
		if [[ -s diff1.txt || -s diff2.txt ]]; then
			echo "IDs in .bim file are not the same as the IDs in the PCA files. Check that you have the right files and run again."
			exit 1
		fi

		~/flashpca/flashpca \
			--bfile "~{output_basename}" \
			--numthreads "~{nthreads}" \
			--project \
			--inmeansd meansd.txt \
			--outproj projections.txt \
			--inload loadings.txt \
			-v
	>>>

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
	}

	output {
		File projections = "~{output_basename}_projections.txt"
	}
}

task MakePCAPlot {
	input {
		File population_pcs
		File target_pcs
		String output_basename
		String docker_image
		Int disk_size = 100
	}

	command <<<
		Rscript -<< "EOF"
			ibrary(dplyr)
			library(readr)
			library(ggplot2)

			population_pcs <- read_tsv("~{population_pcs}")
			target_pcs <- read_tsv("~{target_pcs}")

			ggplot(population_pcs, aes(x=PC1, y=PC2, color="Population")) +
				geom_point(size=0.1, alpha=0.1) +
				geom_point(data = target_pcs, aes(x=PC1, y=PC2, color="Target")) +
				labs(x="PC1", y="PC2") +
				theme_bw()

			ggsave(filename = "~{output_basename}_pca_plot.png", dpi=300, width = 6, height = 6)

		EOF
	>>>

	output {
		File pca_plot = "~{output_basename}_pca_plot.png"
	}

	runtime {
		docker: "~{docker_image}"
 		disks: "local-disk " + disk_size + " HDD"
	}
}

task AdjustScores {
	input {
		File model_parameters
		File pcs
		File scores
		String output_basename
		String docker_image
		Int disk_size = 100
		Int mem_size = 2
	}

	command <<<
		Rscript -<< "EOF"
			library(dplyr)
			library(readr)

			# Read in model parameters
			params_tibble <- read_tsv("~{model_parameters}")
			params <- params_tibble %>% pull(value)

			# Linear transformation to predict variance
			f_sigma2 <- function(t, theta) {
				PC1 = t %>% pull(PC1)
				PC2 = t %>% pull(PC2)
				PC3 = t %>% pull(PC3)
				PC4 = t %>% pull(PC4)
				PC5 = t %>% pull(PC5)
				sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
			}


			# Linear transformation to predict mean
			f_mu <- function(t, theta) {
				PC1 = t %>% pull(PC1)
				PC2 = t %>% pull(PC2)
				PC3 = t %>% pull(PC3)
				PC4 = t %>% pull(PC4)
				PC5 = t %>% pull(PC5)
				mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
			}

			scores = inner_join(read_tsv("~{pcs}"), read_tsv("~{scores}"), by=c("IID" = "#IID"))

			adjusted_scores <- scores %>% mutate(
				adjusted_score = (SCORE1_SUM - f_mu(scores, params[1:5]))/sqrt(f_sigma2(scores, params[6:10]))
				)
			adjusted_scores <- adjusted_scores %>% mutate(percentile=pnorm(adjusted_score,0))

			# Return array scores
			write_tsv(adjusted_scores, "~{output_basename}_adjusted_scores.tsv")
		EOF
	>>>

	output {
		File adjusted_scores = "~{output_basename}_adjusted_scores.tsv"
	}

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
	}
}