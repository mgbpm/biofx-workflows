version 1.0

# The following PRS tasks adapted from:
# https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/ScoringTasks.wdl
# https://raw.githubusercontent.com/broadinstitute/palantir-workflows/refs/heads/main/ImputationPipeline/PCATasks.wdl

# plink chromosome encoding rules: https://www.cog-genomics.org/plink/2.0/data#irreg_output
task DetermineChromosomeEncoding {
	input {
		File condition_zip_file
		String docker_image
		Int disk_size = ceil(size(condition_zip_file, "GB")) + 10
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<

		var_weights_file=$(tar -tf "~{condition_zip_file}" | grep harmonized_weights | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${var_weights_file}
		mv ${var_weights_file} harmonized_weights.txt

		python3 -c '
code = "MT"
with open("harmonzied_weights.txt") as weights_file:
	chroms = {s.split("\t")[0].split(":")[0] for i, s in enumerate(weights_file) if i > 0}
	if any("chr" in c for c in chroms):
		if "chrM" in chroms:
			code = "chrM"
		else:
			code = "chrMT"
	elif "M" in chroms:
		code = "M"
with open("chr_encode_out.txt", "w") as write_code_file:
	write_code_file.write(f"{code}\n")'
	>>>

	runtime {
		docker : "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + "GB"
        preemptible: preemptible
	}

	output {
		String chr_encoding = read_string("chr_encode_out.txt")
	}
}

task ScoreVCF {
	input {
		File input_vcf
		String chromosome_encoding
		File condition_zip_file
		String output_basename
		String docker_image
		Int disk_size =  ceil(size(input_vcf, "GB") * 3) + ceil(size(condition_zip_file, "GB")) + 20
		Int base_mem = 16
		Int mem_size = base_mem + 2
		Int plink_mem = ceil(base_mem * 0.75 * 1000)
		Int preemptible = 1
	}
  
	command <<<
		set -euxo pipefail

		mkdir -p WORK
		sites_file=$(tar -tf "~{condition_zip_file}" | grep .sscore.vars | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${sites_file}
		mv ${sites_file} WORK/sites.sscore.vars

		var_weights_file=$(tar -tf "~{condition_zip_file}" | grep harmonized_weights | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${var_weights_file}
		mv ${var_weights_file} harmonized_weights.txt

		/plink2 --score harmonzied_weights \
			header ignore-dup-ids list-variants no-mean-imputation \
			cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums \
			--set-all-var-ids "@:#:\\$1:\\$2" \
			-vcf "~{input_vcf}" \
			dosage=DS \
			--new-id-max-allele-len 1000 missing \
			--extract WORK/sites.sscore.vars \
			--out "~{output_basename}" \
			--memory "~{plink_mem}" \
			--output-chr "~{chromosome_encoding}"
	>>>

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
		preemptible: preemptible
	}

	output {
		File score = "~{output_basename}.sscore"
		File log = "~{output_basename}.log"
		File sites_scored = "~{output_basename}.sscore.vars"
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
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

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
		preemptible: preemptible
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
		File condition_zip_file
		String output_basename
		String docker_image
		Int disk_size = 400
		Int mem_size = 8
		Int nthreads = 16
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		mkdir -p WORK

		cp "~{bim}" WORK/"~{output_basename}.bim"
		cp "~{bed}" WORK/"~{output_basename}.bed"
		cp "~{fam}" WORK/"~{output_basename}.fam"

		pc_loadings_file=$(tar -tf "~{condition_zip_file}" | grep .pc.loadings | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${pc_loadings_file}
		mv ${pc_loadings_file} WORK/loadings.txt

		meansd_file=$(tar -tf "~{condition_zip_file}" | grep .pc.meansd | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${meansd_file}
		mv ${meansd_file} WORK/meansd.txt

		# Check if .bim file, pc loadings, and pc meansd files have the same IDs
		# 1. extract IDs, removing first column of .bim file and first rows of the pc files
		awk '{print $2}' WORK/"~{output_basename}.bim" > WORK/bim_ids.txt
		awk '{print $1}' WORK/loadings.txt | tail -n +2 > WORK/pcloadings_ids.txt
		awk '{print $1}' WORK/meansd.txt | tail -n +2 > WORK/meansd_ids.txt

		diff WORK/bim_ids.txt WORK/pcloadings_ids.txt > WORK/diff1.txt
		diff WORK/bim_ids.txt WORK/meansd_ids.txt > WORK/diff2.txt
		diff WORK/pcloadings_ids.txt WORK/meansd_ids.txt > WORK/diff3.txt

		if [[ -s WORK/diff3.txt ]];then
			echo "PC loadings file and PC means file do not contain the same IDs. Check your input files and run again."
			exit 1
		fi

		# Check if diff files are not empty
		if [[ -s WORK/diff1.txt || -s WORK/diff2.txt ]]; then
			echo "IDs in .bim file are not the same as the IDs in the PCA files. Check that you have the right files and run again."
			exit 1
		fi

		~/flashpca/flashpca \
			--bfile "~{output_basename}" \
			--numthreads "~{nthreads}" \
			--project \
			--inmeansd WORK/meansd.txt \
			--outproj "~{output_basename}_projections.txt" \
			--inload WORK/loadings.txt \
			-v
	>>>

	runtime {
		docker: "~{docker_image}"
		disks: "local-disk " + disk_size + " HDD"
		memory: mem_size + " GB"
		preemptible: preemptible
	}

	output {
		File projections = "~{output_basename}_projections.txt"
	}
}

task MakePCAPlot {
	input {
		File condition_zip_file
		File target_pcs
		String output_basename
		String docker_image
		Int disk_size = 100
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		pcs_file=$(tar -tf "~{condition_zip_file}" | grep .pc$ | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${pcs_file}
		mv ${pcs_file} population_pcs.txt
		
		Rscript -<< "EOF"
			library(dplyr)
			library(readr)
			library(ggplot2)

			population_pcs <- read_tsv(population_pcs.txt)
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
		memory: mem_size + " GB"
		preemptible: preemptible
	}
}

task AdjustScores {
	input {
		File condition_zip_file
		File pcs
		File scores
		String output_basename
		String docker_image
		Int disk_size = 100
		Int mem_size = 2
		Int preemptible = 1
	}

	command <<<
		set -euxo pipefail

		model_params_file=$(tar -tf "~{condition_zip_file}" | grep fitted_model_params | cut -d "/" -f 2)
		tar -xf "~{condition_zip_file}" ${model_params_file}
		mv ${model_params_file} model_params.tsv

		Rscript -<< "EOF"
			library(dplyr)
			library(readr)

			# Read in model parameters
			params_tibble <- read_tsv(model_params.tsv)
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
		preemptible: preemptible
	}
}