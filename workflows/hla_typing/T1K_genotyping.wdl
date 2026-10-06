version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/steps/Utilities.wdl"

workflow T1kHlaGenotyping {
    input {
        File   input_sample
        File   input_sample_idx
        File   ref_fasta        = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File   ref_fai          = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        String t1k_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/t1k:0.0.1"
    }

    String format = if (basename(input_sample) == basename(input_sample, ".bam") + ".bam") then "bam"
                    else if (basename(input_sample) == basename(input_sample, ".cram") + ".cram") then "cram"
                    else "unsupported"

    if (format == "unsupported") {
        call Utilities.FailTask as UnsupportedFormatError {
            input:
                error_message = "The input sample file is neither a BAM nor CRAM."
        }
    }

    if (format == "bam" || format == "cram") {
        call RunT1kTask {
            input:
                input_format = format,
                input_file   = input_sample,
                input_index  = input_sample_idx,
                ref_fasta    = ref_fasta,
                ref_fai      = ref_fai,
                docker_image = t1k_docker_image
        }
    }

    output {
        File allele_tsv   = select_first([RunT1kTask.allele_tsv])
		File allele_vcf   = select_first([RunT1kTask.allele_vcf])
        File genotype_tsv = select_first([RunT1kTask.genotype_tsv])
    }
}

task RunT1kTask {
    input {
        String input_format
        File   input_file
        File   input_index
        File   ref_fasta
        File   ref_fai
        String preset_type     = "hla"
        String output_basename = basename(input_file, ".~{input_format}")
        String docker_image
        Int    mem_size        = ceil(16 + size(input_file, "GB") * 0.5)
        Int    addldisk        = 8
        Int    preemptible     = 1
    }

    Int sample_size        = ceil(size(input_file, "GB") + size(input_index, "GB"))
    Int ref_size           = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB"))
    Int storage_multiplier = if input_format == "CRAM" then 4 else 2
    Int final_disk_size    = (sample_size * storage_multiplier) + ref_size + addldisk

    command <<<
        set -euxo

        mkdir --parents "OUTPUT"
        mkdir --parents "TEMP"

        if [[ '~{input_format}' == 'CRAM' ]]
        then
            NTHREADS="$( nproc )"

            samtools view --bam              \
                --threads   "${NTHREADS}"    \
                --reference "~{ref_fasta}"   \
                --output    "TEMP/input.bam" \
                "~{input_file}"

            samtools index --bai           \
                --threads "${NTHREADS}"    \
                --output  "TEMP/input.bai" \
                "TEMP/input.bam"

            t1k "~{preset_type}" "TEMP/input.bam" "OUTPUT/~{output_basename}"
        else
            t1k "~{preset_type}" "~{input_file}" "OUTPUT/~{output_basename}"
        fi
    >>>

    output {
        Array[File] output_fa    = glob("~{output_basename}_aligned*.fa")
        Array[File] output_fq    = glob("~{output_basename}_candidate*.fq")
        File        allele_tsv   = "~{output_basename}_allele.tsv"
        File        genotype_tsv = "~{output_basename}_genotype.tsv"
        File        allele_vcf   = "~{output_basename}_allele.vcf"
    }

    runtime {
        docker:      docker_image
        disks:       "local-disk ~{final_disk_size} SSD"
        memory:      "~{mem_size}GB"
        preemptible: preemptible
    }
}