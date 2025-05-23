version 1.0

workflow Glimpse2Imputation {
    input {
        # List of files, one per line
        File reference_chunks

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index
        String output_basename

        File ref_dict
        String af_cutoff
        File gnomadAF_ref_vcf

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Boolean keep_monomorphic_ref_sites = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size

        Boolean collect_qc_metrics = true
        
        Int preemptible = 1
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String docker_extract_num_sites_from_reference_chunk = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        File? monitoring_script
    }

    scatter (reference_chunk in read_lines(reference_chunks)) {
        call GetNumberOfSitesInChunk {
            input:
                reference_chunk = reference_chunk,
                docker = docker_extract_num_sites_from_reference_chunk
        }

        Int n_rare = GetNumberOfSitesInChunk.n_rare
        Int n_common = GetNumberOfSitesInChunk.n_common

        if (defined(input_vcf)) {
            call CountSamples {
                input:
                    vcf = select_first([input_vcf])
            }
        }

        Int nCrams = if defined(crams) then length(select_first([crams])) else 0

        #Int n_samples = select_first([CountSamples.nSamples, length(select_first([crams]))])
        Int n_samples = select_first([CountSamples.nSamples, nCrams])

        call SelectResourceParameters {
            input:
                n_rare = n_rare,
                n_common = n_common,
                #n_samples = select_first([CountSamples.nSamples, length(select_first([crams]))])
                n_samples = n_samples
        }

        if (SelectResourceParameters.memory_gb > 256 || SelectResourceParameters.request_n_cpus > 32) {
            # force failure if we're accidently going to request too much resources and spend too much money
            Int safety_check_memory_gb = -1
            Int safety_check_n_cpu = -1
        }
        call GlimpsePhase {
            input:
                reference_chunk = reference_chunk,
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                impute_reference_only_variants = impute_reference_only_variants,
                keep_monomorphic_ref_sites = keep_monomorphic_ref_sites,
                n_burnin = n_burnin,
                n_main = n_main,
                effective_population_size = effective_population_size,
                call_indels = call_indels,
                crams = crams,
                cram_indices = cram_indices,
                sample_ids = sample_ids,
                fasta = fasta,
                fasta_index = fasta_index,
                preemptible = preemptible,
                docker = docker,
                cpu = select_first([safety_check_n_cpu, SelectResourceParameters.request_n_cpus]),
                mem_gb = select_first([safety_check_memory_gb, SelectResourceParameters.memory_gb]),
                monitoring_script = monitoring_script
        }
    }

    call GlimpseLigate {
        input:
            imputed_chunks = GlimpsePhase.imputed_vcf,
            imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
            output_basename = output_basename,
            ref_dict = ref_dict,
            preemptible = preemptible,
            docker = docker,
            monitoring_script = monitoring_script
    }

    call Filter_af {
        input:
            input_vcf = GlimpseLigate.imputed_vcf,
            gnomadAF_ref_vcf = gnomadAF_ref_vcf,
            af_cutoff = af_cutoff,
            output_basename = output_basename,
    }

    if (collect_qc_metrics) {
        call CollectQCMetrics {
            input:
                imputed_vcf = GlimpseLigate.imputed_vcf,
                output_basename = output_basename,
                monitoring_script = monitoring_script
        }
    }

    output {
        File imputed_afFiltered_vcf = Filter_af.imputed_afFiltered_vcf
        File imputed_afFiltered_vcf_index = Filter_af.imputed_afFiltered_vcf_index
        File? qc_metrics = CollectQCMetrics.qc_metrics

        Array[File?] glimpse_phase_monitoring = GlimpsePhase.monitoring
        File? glimpse_ligate_monitoring = GlimpseLigate.monitoring
    }
}

task GlimpsePhase {
    input {
        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index
        File reference_chunk

        Boolean impute_reference_only_variants
        Boolean keep_monomorphic_ref_sites
        Boolean call_indels
        Int? n_burnin
        Int? n_main
        Int? effective_population_size

        Int mem_gb = 16
        Int cpu = 8
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + size(reference_chunk, "GiB") + 10)
        Int preemptible = 1
        Int max_retries = 2
        String docker
        File? monitoring_script
    }

    parameter_meta {
        crams: {
                        localization_optional: true
                    }
        cram_indices: {
                        localization_optional: true
                    }
    }

    String bam_file_list_input = if defined(crams) then "--bam-list crams.list" else ""
    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cram_paths=( ~{sep=" " crams} )
        sample_ids=( ~{sep=" " sample_ids} )

        duplicate_cram_filenames=$(printf "%s\n" "${cram_paths[@]}" | xargs -I {} basename {} | uniq -d)
        if [ ! -z "$duplicate_cram_filenames" ]; then
            echo "ERROR: The input CRAMs contain multiple files with the same basename, which leads to an error due to the way that htslib is implemented. Duplicate filenames:"
            printf "%s\n" "${duplicate_cram_filenames[@]}"
            exit 1
        fi

        for i in "${!cram_paths[@]}" ; do
            echo -e "${cram_paths[$i]} ${sample_ids[$i]}" >> crams.list
        done

        cmd="/bin/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
        ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} ~{if call_indels then "--call-indels" else ""} \
        ~{if keep_monomorphic_ref_sites then "--keep-monomorphic-ref-sites" else ""} \
        ~{"--burnin " + n_burnin} ~{"--main " + n_main} \
        ~{"--ne " + effective_population_size} \
        ~{bam_file_list_input} \
        ~{"--fasta " + fasta} \
        --checkpoint-file-out checkpoint.bin"

        if [ -s "checkpoint.bin" ]; then
            cmd="$cmd --checkpoint-file-in checkpoint.bin" 
        fi

        eval $cmd
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        checkpointFile: "checkpoint.bin"
    }

    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
        File? monitoring = "monitoring.log"
    }
}

task GlimpseLigate {
    input {
        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices
        String output_basename
        File ref_dict

        Int mem_gb = 24
        Int cpu = 8
        Int disk_size_gb = ceil(2.2 * size(imputed_chunks, "GiB") + 100)
        Int preemptible = 1
        Int max_retries = 1
        String docker
        File? monitoring_script
    }

    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."
        
        /bin/GLIMPSE2_ligate --input ~{write_lines(imputed_chunks)} --output ligated.vcf.gz --threads ${NPROC}

        # sort ligated.vcf first otherwise it does not run other bcftools commands
        # change ID and deduplicate it
        bcftools sort ligated.vcf.gz -Ou | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou | bcftools norm -d both -Oz -o ligated_cleaned.vcf.gz

        # Set correct reference dictionary
        bcftools view -h --no-version ligated.vcf.gz > old_header.vcf
        java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD ~{ref_dict} -O new_header.vcf        
        bcftools reheader -h new_header.vcf -o ~{output_basename}.imputed.vcf.gz ligated_cleaned.vcf.gz
        tabix ~{output_basename}.imputed.vcf.gz

        # Remove intermediary files. Note: Sometimes, ligated.vcf.gz.csi cannot be created due to "[E::hts_idx_push] Chromosome blocks not continuous"
        rm ligated.vcf.gz
        rm ligated_cleaned.vcf.gz
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
        File? monitoring = "monitoring.log"
    }
}

task CollectQCMetrics {
    input {
        File imputed_vcf
        String output_basename
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    parameter_meta {
        imputed_vcf: {
                        localization_optional: true
                    }
    }

    Int disk_size_gb = 100
    
    command <<<
        set -euo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cat <<'EOF' > script.py
import hail as hl
import pandas as pd

# Calculate metrics
hl.init(default_reference='GRCh38', idempotent=True)
vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
qc = hl.sample_qc(vcf)
qc.cols().flatten().rename({'sample_qc.' + col: col for col in list(qc['sample_qc'])}).export('~{output_basename}.qc_metrics.tsv')
EOF
        python3 script.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File qc_metrics = "~{output_basename}.qc_metrics.tsv"
        File? monitoring = "monitoring.log"
    }
}

task GetNumberOfSitesInChunk {
    input {
        File reference_chunk

        String docker
        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(size(reference_chunk, "GiB") + 10)
        Int preemptible = 1
        Int max_retries = 1
    }

    command <<<
        set -xeuo pipefail
        /bin/GLIMPSE2_extract_num_sites_from_reference_chunk ~{reference_chunk} > n_sites.txt
        cat n_sites.txt
        grep "Lrare" n_sites.txt | sed 's/Lrare=//' > n_rare.txt
        grep "Lcommon" n_sites.txt | sed 's/Lcommon=//' > n_common.txt
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        Int n_rare = read_int("n_rare.txt")
        Int n_common = read_int("n_common.txt")
    }
}

task CountSamples {
  input {
    File vcf

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = 10 + ceil(size(vcf, "GiB"))
  }

  command <<<
    bcftools query -l ~{vcf} | wc -l
  >>>

  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
  }
  output {
    Int nSamples = read_int(stdout())
  }
}

task SelectResourceParameters {
    input {
        Int n_rare
        Int n_common
        Int n_samples
    }

    command <<<
        python3 << EOF
        import math
        n_rare = ~{n_rare}
        n_common = ~{n_common}
        n_samples = ~{n_samples}
        n_sites = n_common + n_rare

        # try to keep expected runtime under 4 hours, but don't ask for more than 32 cpus, or 256 GB memory
        estimated_needed_threads = min(math.ceil(5e-6*n_sites*n_samples/240), 32)
        estimated_needed_memory_gb = min(math.ceil((800e-3 + 0.97e-6 * n_rare * estimated_needed_threads + 14.6e-6 * n_common * estimated_needed_threads + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples))), 256)
        #estimated_needed_memory_gb = min(math.ceil((1 + 1.5e-6 * n_rare * estimated_needed_threads + 25e-6 * n_common * estimated_needed_threads + 10e-9 * (n_rare + n_common) * n_samples + 25e-3 * n_samples + 3e-6*(n_rare + n_common)*math.log(n_samples))), 256)

        # recalc allowable threads, may be some additional threads available due to rounding memory up
        threads_to_use = max(math.floor((estimated_needed_memory_gb - (800e-3 + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples)))/(0.97e-6 * n_rare + 14.6e-6 * n_common)), 1) 
        estimated_needed_memory_gb = math.ceil(2 * estimated_needed_memory_gb)

        with open("n_cpus_request.txt", "w") as f_cpus_request:
            f_cpus_request.write(f'{int(threads_to_use)}')

        with open("memory_gb.txt", "w") as f_mem:
            f_mem.write(f'{int(estimated_needed_memory_gb)}')
        EOF
    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
    }

    output {
        Int memory_gb = read_int("memory_gb.txt")
        Int request_n_cpus = read_int("n_cpus_request.txt")
    }
}

task Filter_af {
    input{
        String input_vcf
        String af_cutoff = ">=0.0001"
        String gnomadAF_ref_vcf
        String output_basename

        Int mem_gb = 16
        Int cpu = 4
        Int disk_size_gb = ceil(2.5 * (size(input_vcf, "GiB") + size(gnomadAF_ref_vcf, "GiB"))) + 10
        Int preemptible = 1
        Int max_retries = 1
        String bcftools_docker = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
    }

    String input_vcf_name = basename(input_vcf)
    String input_vcf_index = basename(input_vcf) + ".tbi"
    #String input_vcf_filtered_name = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + "_filtered.vcf.gz"
    #String input_vcf_filtered_name_index = input_vcf_filtered_name + ".tbi"

    String gnomadAF_ref_vcf_name = basename(gnomadAF_ref_vcf)
    String gnomadAF_ref_vcf_index = basename(gnomadAF_ref_vcf) + ".tbi"
    String gnomadAF_ref_vcf_filtered_name = sub(basename(gnomadAF_ref_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + "_filtered.vcf.gz"

    File input_vcf_file = input_vcf
    File input_vcf_file_index = input_vcf_file + ".tbi"
    File gnomadAF_ref_vcf_file = gnomadAF_ref_vcf
    File gnomadAF_ref_vcf_file_index = gnomadAF_ref_vcf_file + ".tbi"

    command <<<
        set -euo pipefail

        mkdir vcf_dir
        mkdir sort_tmp_dir
        ln -s ~{input_vcf_file} vcf_dir/~{input_vcf_name}
        ln -s ~{input_vcf_file_index} vcf_dir/~{input_vcf_index}
        ln -s ~{gnomadAF_ref_vcf_file} vcf_dir/~{gnomadAF_ref_vcf_name}
        ln -s ~{gnomadAF_ref_vcf_file_index} vcf_dir/~{gnomadAF_ref_vcf_index}

        bcftools filter -i "INFO/VEP_gnomad4.1_joint_frequency_AF_grpmax_joint ~{af_cutoff}" \
                 vcf_dir/~{gnomadAF_ref_vcf_name} \
                 -Oz \
                 -o ~{gnomadAF_ref_vcf_filtered_name}
        bcftools index -ft ~{gnomadAF_ref_vcf_filtered_name}

        bcftools isec -p '.' -w1 -Oz vcf_dir/~{input_vcf_name} ~{gnomadAF_ref_vcf_filtered_name}

        #mv 0002.vcf.gz ~{output_basename}.imputed.filtered.vcf.gz
        #mv 0002.vcf.gz.tbi ~{output_basename}.imputed.filtered.vcf.gz.tbi

        bcftools sort --temp-dir sort_tmp_dir -Oz -o ~{output_basename}.imputed.filtered.vcf.gz 0002.vcf.gz
        bcftools index -ft ~{output_basename}.imputed.filtered.vcf.gz
    >>>

    runtime {
        docker: bcftools_docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_afFiltered_vcf = "~{output_basename}.imputed.filtered.vcf.gz"
        File imputed_afFiltered_vcf_index = "~{output_basename}.imputed.filtered.vcf.gz.tbi"
    }
}
