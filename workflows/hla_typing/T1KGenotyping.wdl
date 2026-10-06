version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/steps/Utilities.wdl"

workflow T1kHlaGenotyping {
    input {
        File    input_sample
        File    input_sample_idx
        File?   hla_intervals
        String? google_project_id
        File    ref_fasta        = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File    ref_fai          = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File    ref_dict         = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
        String  t1k_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/t1k:0.0.1"
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
        if (defined(hla_intervals)) {
            call FilterBamToHLA {
                input:
                    original_bam     = input_sample,
                    original_bam_idx = input_sample_idx,
                    ref_fasta        = ref_fasta,
                    ref_fai          = ref_fai,
                    ref_dict         = ref_dict,
                    hla_intervals    = select_first([hla_intervals]),
                    google_project   = google_project_id
            }

            String new_format = "bam"
        }
        
        call RunT1kTask {
            input:
                input_format = select_first([new_format, format]),
                input_file   = select_first([FilterBamToHLA.hla_bam, input_sample]),
                input_index  = select_first([FilterBamToHLA.hla_bam_idx, input_sample_idx]),
                ref_fasta    = ref_fasta,
                ref_fai      = ref_fai,
                docker_image = t1k_docker_image
        }
    }

    output {
        File? hla_bam = FilterBamToHLA.hla_bam
        File? hla_bai = FilterBamToHLA.hla_bam_idx
        File  allele_tsv   = select_first([RunT1kTask.allele_tsv])
		File  allele_vcf   = select_first([RunT1kTask.allele_vcf])
        File  genotype_tsv = select_first([RunT1kTask.genotype_tsv])
    }
}

task FilterBamToHLA {
    input {
        File    original_bam       # this can be a BAM or CRAM
        File    original_bam_idx
        File    ref_fasta          # GATK PrintReads requires a reference for CRAMs
        File    ref_fai
        File    ref_dict
        File    hla_intervals
        String? google_project
        Int     cpu          = 2
        Int     num_threads  = 4
        Int     mem_size     = 4
        Int     addldisk     = 100
        Int     boot_disk_gb = 10
        Int     preemptible  = 1
        String  docker_image = "us.gcr.io/broad-gatk/gatk"
    }

    Int ref_size    = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB"))
    Int sample_size = ceil(size(original_bam, "GB") + size(original_bam_idx, "GB"))
    Int disk_size   = addldisk + ref_size + ceil(size(original_bam, "GB"))

    parameter_meta{
        hla_intervals: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fai: {localization_optional: true}
        ref_dict: {localization_optional: true}
        original_bam: {localization_optional: true}
        original_bam_idx: {localization_optional: true}
    }

    command <<<
        gatk PrintReads -R ~{ref_fasta} -I ~{original_bam} -L ~{hla_intervals} -O hla-unsorted.bam \
            ~{if select_first([google_project, ""]) != "" then "--gcs-project-for-requester-pays " + select_first([google_project, ""]) else ""}

        gatk ValidateSamFile -I hla-unsorted.bam

        samtools sort -@ ~{num_threads} hla-unsorted.bam > hla.bam

        gatk BuildBamIndex -I hla.bam
    >>>

    runtime {
        docker:         docker_image
        bootDiskSizeGb: boot_disk_gb
        memory:         mem_size + " GB"
        disks:          "local-disk " + disk_size + " SSD"
        preemptible:    preemptible
        cpu:            cpu
    }

    output {
        File hla_bam     = "hla.bam"
        File hla_bam_idx = "hla.bai"
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
        String output_basename = sub(basename(input_file), "\\.(cram|CRAM|bam|BAM)$", "")
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

    runtime {
        docker:      docker_image
        disks:       "local-disk ~{final_disk_size} SSD"
        memory:      "~{mem_size}GB"
        preemptible: preemptible
    }

    output {
        Array[File] output_fa    = glob("~{output_basename}_aligned*.fa")
        Array[File] output_fq    = glob("~{output_basename}_candidate*.fq")
        File        allele_tsv   = "~{output_basename}_allele.tsv"
        File        genotype_tsv = "~{output_basename}_genotype.tsv"
        File        allele_vcf   = "~{output_basename}_allele.vcf"
    }
}