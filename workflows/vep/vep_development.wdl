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
        File? custom_database
        # VEP command string with custom database load instructions
        String? custom_datatbase_config
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
        }
    }

    output {
        File? database_output = VEPDatabaseTask.output_vcf_file
        File? cache_output = VEPCacheTask.output_vcf_file
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
    }

    Int disk_size = ceil(size(input_vcf, "GB") + size(cache_file, "GB") * 2 ) + 20 + select_first([extra_disk_gb, 0])
    Int machine_mem_gb = select_first([mem_gb, 20])
    Int cpu_count = select_first([machine_cpu, 10])
    Int thread_count = cpu_count * 2

    command <<<
        set -euxo pipefail

        # Assess allocated disk and memory resources
        echo "Current memory requested: ~{machine_mem_gb} GB"
        echo "Current disk requested: ~{disk_size} GB"
        echo "Current cpu count: ~{cpu_count}"
        echo "Current thread count: ~{thread_count}"

        # Ensure the destination directory exists
        mkdir -p /cromwell_root/.vep

        #investigate available disks
        df -h

        # Unpack the tar.gz file into /cromwell_root/.vep
        tar xzf ~{cache_file} -C /cromwell_root/.vep

        #investigate available disks after unpacking
        df -h

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
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD" 
        memory: machine_mem_gb + " GB"
        cpus: ~{cpu_count}
    }

    output {
        File output_vcf_file = "~{output_name}.vcf.gz"
    }
}