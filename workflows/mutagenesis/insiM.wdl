version 1.0

workflow insiMWorkflow {
    input {
        File input_bam
        File input_bai
        File genome_ref_fasta
        File genome_ref_amb
        File genome_ref_ann
        File genome_ref_bwt
        File genome_ref_pac
        File genom_ref_sa
        String assay_type
        String mutation_type
        String output_basename
        String vaf = ""
        String inslen = ""
        String insseq = ""
        String? amp_read_length
        String docker_image
        Array[MutationTargets] targets
        Array[AmpliconBed] ampbed = []
    }

    if (assay_type == 'amplicon') {
        call MutateAmpliconAssayTask as MutateAmplicon {
            input:
                input_bam = input_bam,
                input_bai = input_bai,
                assay_type = assay_type,
                targets = targets,
                mutation_type = mutation_type,
                genome_ref_fasta = genome_ref_fasta,
                amp_read_length = amp_read_length,
                vaf = vaf,
                inslen = inslen,
                insseq = insseq,
                ampbed = ampbed,
                output_basename = output_basename,
                docker_image = docker_image
        }

        call ConvertToSamTask as AmpliconSam {
            input:
                genome_ref_fasta = genome_ref_fasta,
                genome_ref_amb = genome_ref_amb,
                genome_ref_ann = genome_ref_ann,
                genome_ref_bwt = genome_ref_bwt,
                genome_ref_pac = genome_ref_pac,
                genom_ref_sa = genom_ref_sa,
                in_fastq1 = MutateAmplicon.fastq1,
                in_fastq2 = MutateAmplicon.fastq2,
                input_bam = input_bam,
                output_basename = output_basename
        }

        call ConvertToBamTask as AmpliconBam {
            input:
                output_basename = output_basename,
                mut_sam = AmpliconSam.mut_sam,
                input_bam = input_bam
        }

        call IndexBamTask as AmpliconIndexBam {
            input:
                output_basename = output_basename,
                mut_dedup_bam = AmpliconBam.mut_dedup_bam,
                input_bam = input_bam
        }
    }

    if (assay_type == 'capture') {
        call MutateCaptureAssayTask as MutateCapture {
            input:
                input_bam = input_bam,
                input_bai = input_bai,
                assay_type = assay_type,
                targets = targets,
                mutation_type = mutation_type,
                genome_ref_fasta = genome_ref_fasta,
                vaf = vaf,
                inslen = inslen,
                insseq = insseq,
                output_basename = output_basename,
                docker_image = docker_image
        }

        call ConvertToSamTask as CaptureSam {
            input:
                genome_ref_fasta = genome_ref_fasta,
                genome_ref_amb = genome_ref_amb,
                genome_ref_ann = genome_ref_ann,
                genome_ref_bwt = genome_ref_bwt,
                genome_ref_pac = genome_ref_pac,
                genom_ref_sa = genom_ref_sa,
                in_fastq1 = MutateCapture.fastq1,
                in_fastq2 = MutateCapture.fastq2,
                input_bam = input_bam,
                output_basename = output_basename
        }

        call ConvertToBamTask as CaptureBam {
            input:
                output_basename = output_basename,
                mut_sam = CaptureSam.mut_sam,
                input_bam = input_bam
        }

        call IndexBamTask as CaptureIndexBam {
            input:
                output_basename = output_basename,
                mut_dedup_bam = CaptureBam.mut_dedup_bam,
                input_bam = input_bam
        }
    }

    output {
        File fastq1 = select_first([MutateCapture.fastq1, MutateAmplicon.fastq1])
        File fastq2 = select_first([MutateCapture.fastq2, MutateAmplicon.fastq2])
        File vcf = select_first([MutateCapture.vcf, MutateAmplicon.vcf])
        File mut_sam = select_first([CaptureSam.mut_sam, AmpliconSam.mut_sam])
        File mut_dedup_bam = select_first([CaptureBam.mut_dedup_bam, AmpliconBam.mut_dedup_bam])
        File mut_sorted_bam = select_first([CaptureBam.mut_sorted_bam, AmpliconBam.mut_sorted_bam])
        File metrics_file = select_first([CaptureBam.metrics_file, AmpliconBam.metrics_file])
        File mut_dedup_bai = select_first([CaptureIndexBam.mut_dedup_bai, AmpliconIndexBam.mut_dedup_bai])
    }
}

struct MutationTargets {
    String chr
    Int start
    Int end
    String? type
    Float? vaf
    String? ins_seq
    Int? ins_len
}

struct AmpliconBed {
    String chr
    Int start
    Int end
}

# add mutations to amplicon assay
task MutateAmpliconAssayTask {
    input {
        File input_bam
        File input_bai
        String assay_type
        Array[MutationTargets] targets
        Array[AmpliconBed] ampbed = []
        String mutation_type
        File genome_ref_fasta
        String vaf = ""
        String inslen = ""
        String insseq = ""
        String? amp_read_length
        String output_basename
        String docker_image
    }

    command <<<
        set -euxo pipefail
        $MGBPMBIOFXPATH/insiM/create_insiM_input.py \
            --assay "~{assay_type}" \
            --input-bam "~{input_bam}" \
            --target-json-file "~{write_json(targets)}" \
            --mutation-type "~{mutation_type}" \
            --genome-ref "~{genome_ref_fasta}" \
            --config-file "config_amplicon.txt" \
            --ampliconbed-json-file "~{write_json(ampbed)}" \
            --read-length "~{amp_read_length}" \
            --vaf "~{vaf}" \
            --inslen "~{inslen}" \
            --insseq "~{insseq}" \
            --output-file "~{output_basename}"

        $MGBPMBIOFXPATH/insiM/insiM_v2.0.py -config "config_amplicon.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: "14GB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2.5) + 20) + " HDD"
    }

    output {
        File fastq1 = output_basename + "_R1_001.fastq"
        File fastq2 = output_basename + "_R2_001.fastq"
        File vcf = output_basename + ".insiM.vcf"
    }
}

# add mutation to hybrid capture assay
task MutateCaptureAssayTask {
    input {
        File input_bam
        File input_bai
        String assay_type
        Array[MutationTargets] targets
        String mutation_type
        File genome_ref_fasta
        String vaf = ""
        String inslen = ""
        String insseq = ""
        String output_basename
        String docker_image
    }

    command <<<
        set -euxo pipefail
        $MGBPMBIOFXPATH/insiM/create_insiM_input.py \
            --assay "~{assay_type}" \
            --input-bam "~{input_bam}" \
            --target-json-file "~{write_json(targets)}" \
            --mutation-type "~{mutation_type}" \
            --genome-ref "~{genome_ref_fasta}" \
            --config-file "config_capture.txt" \
            --vaf "~{vaf}" \
            --inslen "~{inslen}" \
            --insseq "~{insseq}" \
            --output-file "~{output_basename}"

        $MGBPMBIOFXPATH/insiM/insiM_v2.0.py -config "config_capture.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: "14GB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2.5) + 20) + " HDD"
    }

    output {
        File fastq1 = output_basename + "_R1_001.fastq"
        File fastq2 = output_basename + "_R2_001.fastq"
        File vcf = output_basename + ".insiM.vcf"
    }
}

# convert FASTQ output from mutation tasks to BAM
task ConvertToSamTask {
    input {
        File genome_ref_fasta
        File genome_ref_amb
        File genome_ref_ann
        File genome_ref_bwt
        File genome_ref_pac
        File genom_ref_sa
        File in_fastq1
        File in_fastq2
        File input_bam
        String output_basename
    }

    command <<<
        set -euxo pipefail

        # align the FASTQ's to reference and produce SAM output
        bwa mem -M "~{genome_ref_fasta}" "~{in_fastq1}" "~{in_fastq2}" \
            > "~{output_basename}"_aligned_reads.sam
    >>>

    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
        memory: "14GB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2.5) + 20) + " HDD"
    }

    output {
        File mut_sam = output_basename + "_aligned_reads.sam"
    }
}

# convert SAM to BAM
task ConvertToBamTask{
    input {
        File mut_sam
        File input_bam
        String output_basename
    }

    command <<<
        set -euxo pipefail

        java -jar /usr/picard/picard.jar SortSam \
            --INPUT "~{mut_sam}" \
            --OUTPUT "~{output_basename}"_sorted_reads.bam \
            --SORT_ORDER coordinate

        java -jar /usr/picard/picard.jar MarkDuplicates \
            --INPUT "~{output_basename}"_sorted_reads.bam \
            --OUTPUT "~{output_basename}"_dedup_reads.bam \
            --METRICS_FILE "~{output_basename}"_metrics.txt
    >>>

    runtime {
        docker: "biocontainers/picard:v1.139_cv3"
        memory: "14GB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2.5) + 20) + " HDD"
    }

    output {
        File mut_sorted_bam = output_basename + "_sorted_reads.bam"
        File mut_dedup_bam = output_basename + "_dedup_reads.bam"
        File metrics_file = output_basename + "_metrics.txt"
    }
}

# index the output BAM
task IndexBamTask {
    input {
        File input_bam
        File mut_dedup_bam
        String output_basename
    }

    command <<<
        set -euxo pipefail

        # index the output mutated bam
        samtools index -b "~{mut_dedup_bam}" "~{output_basename}"_dedup_reads.bam.bai
    >>>

    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2.5) + 20) + " HDD"
    }

    output {
        File mut_dedup_bai = output_basename + "_dedup_reads.bam.bai"
    }
}