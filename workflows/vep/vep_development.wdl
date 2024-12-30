version 1.0

workflow VEPWorkflow {
    input {
        # VCF to annotate with VEP
        File input_vcf
        # Output file name
        String output_name
        # Full path to directory with cache
        File cache_file
        # Cache version number (ex: 111 or 113)
        String cache_version
        # VEP Docker image
        String vep_docker_image = "ensemblorg/ensembl-vep:latest"
        # Choose the annotation source
        Boolean run_database = false
        Boolean run_cache = true
    
        # Custom database file for VEP annotation
        #File? custom_database
        # VEP command string with custom database load instructions
        #String? custom_database_config
    }

    # If a local cache isn't available, run with the online database
    if (run_database) {
        call VEPDatabaseTask {
            input:
                input_vcf = input_vcf,
                output_name = output_name,
                docker_image = vep_docker_image
        }
    }

    # If a local cache is available, run with the cache
    if (run_cache) {


        call VEPCacheTask {
            input:
                input_vcf = input_vcf,
                cache_file = cache_file,
                cache_version = cache_version,
                output_name = output_name,
                docker_image = vep_docker_image
                #custom_database = custom_database
                #custom_database_config = custom_database_config
        }
    }

    output {
        File? database_output = VEPDatabaseTask.output_vcf_file
        File? cache_output_vcf = VEPCacheTask.output_vcf_file
        File? cache_output_vcf = VEPCacheTask.output_vcf_file
    }
    meta {
        allowNestedInputs: true
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
            --output_file "~{output_name}.vcf.gz" --vcf --no_stats \
            --compress_output bgzip \
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
            --assembly GRCh38 \
            --vcf
    >>>

    runtime {
        docker: "~{docker_image}"
    }

    output {
        File output_vcf_file = "~{output_name}.vcf.gz"
    }
}

task VEPCacheTask {
    input {
        File input_vcf
        File cache_file
        String cache_version
        String output_name
        String docker_image
        # Specify an amount of additional disk space to add to the VM
        Int? extra_disk_gb
        # Specify the amount of RAM the VM uses
        Int? mem_gb
        # Number of cpus to use while annotating, default is 10
        Int? machine_cpus

        File? custom_database
        File? custom_database_index
        String? custom_database_config
    }

    Int disk_size = ceil(size(input_vcf, "GB") + size(cache_file, "GB") * 2 ) + 20 + select_first([extra_disk_gb, 0])
    Int machine_mem_gb = select_first([mem_gb, 20])
    Int cpu_count = select_first([machine_cpus, 10])
    Int thread_count = cpu_count * 2

    command <<<
        set -euxo pipefail

        #calculate current cpu count
        real_cpus=$(nproc)

        # Assess allocated disk and memory resources
        echo "Current memory requested: ~{machine_mem_gb} GB"
        echo "Current disk requested: ~{disk_size} GB"
        echo "Current cpu requested: ~{cpu_count}"
        echo "Current thread requested: ~{thread_count}"
        echo "Current nproc count: $real_cpus"

        # Ensure the destination directory exists
        mkdir -p /cromwell_root/.vep

        #investigate available disks
        df -h

        # Unpack the tar.gz file into /cromwell_root/.vep
        tar xzf ~{cache_file} -C /cromwell_root/.vep

        #investigate available disks after unpacking
        df -h

        if [[ -n "~{custom_database}" && -n "~{custom_database_config}" ]]; then
            echo "Both custom database file and custom database config defined."
            use_custom_db="true"
        elif [[ -n "~{custom_database}" ]]; then
            echo "Custom database file provided, missing custom database config."
            use_custom_db="false"
        elif [[ -n "~{custom_database_config}" ]]; then
            echo "Custom database config provided, missing custom database file."
            use_custom_db="false"
        else
            echo "Neither custom database file or custom database config provided."
            use_custom_db="false"
        fi

        if [ ${use_custom_db} == "true" ] ; then
            echo "Localizing additional custom database: ~{custom_database}"
            echo "Custom database configuration: ~{custom_database_config}"
            if [ -n "~{custom_database_index}" ] ; then
                echo "Custom database index also detected: ~{custom_database_index}"
            else
                echo "Custom database index file missing, Exiting"
                exit
            fi

            #run VEP command
            /opt/vep/src/ensembl-vep/vep \
                --cache --dir_cache /cromwell_root/.vep \
                --cache_version "~{cache_version}" \
                --input_file "~{input_vcf}" \
                --output_file "~{output_name}.vcf.gz" \
                --vcf --no_stats \
                --compress_output bgzip \
                --verbose \
                --fork "~{thread_count}" \
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
                --assembly GRCh38 \
                --custom file="~{custom_database}",~{custom_database_config}
        else
            echo "Custom database file or configuration not detected, running normally."
            /opt/vep/src/ensembl-vep/vep \
                --cache --dir_cache /cromwell_root/.vep \
                --cache_version "~{cache_version}" \
                --input_file "~{input_vcf}" \
                --output_file "~{output_name}.vcf.gz" \
                --vcf --no_stats \
                --compress_output bgzip \
                --verbose \
                --fork "~{thread_count}" \
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
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD" 
        memory: machine_mem_gb + " GB"
        cpu: "~{cpu_count}"
    }

    output {
        File output_vcf_file = "~{output_name}.vcf.gz"
    }
}