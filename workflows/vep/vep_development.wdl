version 1.0

workflow VEPWorkflow {
    input {
		# VCF to annotate with VEP
        File input_vcf
		# Output file name
		String output_name
		# Full path to directory with cache
		String cache_location
		# Cache version number (ex: 111 or 113)
		String cache_version
		# VEP Docker image
        String vep_docker_image = "ensemblorg/ensembl-vep:latest"
		# Choose the annotation source
		Boolean run_database = false
		Boolean run_cache = true
    }

	# If a local cache isn't available, run with the online database
	if (run_database) {
    	call VEPDatabaseTask {
    	    input:
    	        input_vcf = input_vcf,
				output_name - output_name,
    	        docker_image = vep_docker_image
    	}
	}

	# If a local cache is available, run with the cache
	if (run_cache) {
    	call VEPCacheTask {
    	    input:
    	        input_vcf = input_vcf,
    	        cache_location = cache_location,
				cache_version = cache_version,
				output_name = output_name,
    	        docker_image = vep_docker_image
    	}
	}

    output {
        File? database_output = VEPDatabaseTask.output_tab_file
        File? cache_output = VEPCacheTask.output_tab_file
    }
}

task VEPDatabaseTask {
    input {
        File input_vcf
		String output_name
        String docker_image
    }

    command <<<
        set -euxo pipefail

        /opt/vep/src/ensembl-vep/vep --database \
            --input_file "~{input_vcf}" \
            --output_file "~{output_name}.txt" --tab --no_stats \
            --verbose \
            --show_ref_allele \
            --symbol \
            --hgvs \
            --canonical \
            --numbers \
            --individual all \
            --check_svs \
            --clin_sig_allele 1 \
            --mane \
            --species homo_sapiens \
			--assembly GRCh38
    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        File output_tab_file = "~{output_name}.txt"
    }
}

task VEPCacheTask {
    input {
        File input_vcf
        String cache_location
		String cache_version
		String output_name
        String docker_image
    }

    command <<<
        set -euxo pipefail

        /opt/vep/src/ensembl-vep/vep \
            --cache --dir_cache "~{cache_location}" --cache_version "~{cache_version}" \
            --input_file "~{input_vcf}" \
            --output_file "~{output_name}.txt" --tab --no_stats \
            --verbose \
            --show_ref_allele \
            --symbol \
            --hgvs \
            --canonical \
            --numbers \
            --individual all \
            --af_gnomade \
            --af_gnomadg \
            --max_af \
            --clin_sig_allele 1 \
            --pubmed \
            --mane \
            --sift b \
            --polyphen b \
            --ccds \
            --domains \
            --regulatory \
            --protein \
            --biotype \
            --af \
            --af_1kg \
            --uniprot \
            --tsl \
            --appris \
            --variant_class \
            --gene_phenotype \
            --mirna \
            --species homo_sapiens \
			--assembly GRCh38
    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        File output_tab_file = "~{output_name}.txt"
    }
}