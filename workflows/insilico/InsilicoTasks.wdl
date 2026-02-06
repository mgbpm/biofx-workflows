version 1.0

struct MutationBED {
    String  chr
    Int     start
    Int     end
    Int?    num_dupes
    Int?    trans_start
    Float?  vaf
    String? mut_type
    String? base
    String? insert_seq
    String? target_site_len
    String? cut_site_motif
    String? trans_chr
    String? trans_end
    String? trans_on_chr
    String? trans_from_chr
}

task CramToBamTask {
    input {
        File   input_cram
        File   input_crai
        File   ref_fasta
        File   ref_fai
        String output_basename = sub(basename(input_cram), ".cram", "")
        String docker_image
        Int    addldisk        = 100 
        Int    mem_size        = 16
        Int    preemptible     = 1
    }

    Int cram_size       = ceil(size(input_cram, "GB") * 3) + ceil(size(input_crai, "GB"))
    Int ref_size        = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB"))
    Int final_disk_size = cram_size + ref_size + addldisk

    command <<<
        set -euxo

        samtools view -b -T "~{ref_fasta}" "~{input_cram}" > "~{output_basename}.bam"
        samtools index -b "~{output_basename}.bam" "~{output_basename}.bam.bai"
    >>>

    runtime {
        docker: docker_image
        disks: "local-disk ~{final_disk_size} SSD"
        memory: "~{mem_size}GB"
        preemptible: preemptible
    }

    output {
        File output_bam = "~{output_basename}.bam"
        File output_bai = "~{output_basename}.bam.bai"
    }
}

task RunBamsurgeonTask {
    input {
        Array[MutationBED] mutation_bed
        File    input_bam
        File    input_bai
        File    ref_fasta
        File    ref_amb
        File    ref_ann
        File    ref_bwt
        File    ref_fai
        File    ref_pac
        File    ref_sa
        File?   cnv_file
        File?   avoid_reads
        File?   donor_bam
        File?   donor_bai
        File?   insert_library
        Int     nums_snvs            = 0
        Int     haplo_size           = 0
        Int     min_depth            = 10
        Int     max_depth            = 2000
        Int     min_mut_reads        = 3
        Int     max_open_files       = 1000
        Int     max_lib_size         = 600
        Int     kmer_size            = 31
        Int     min_contig_gen       = 4000
        Int     mean_insert_size     = 300
        Int     insert_size_stdev    = 70
        Float   cover_diff           = 0.9
        Float   snv_frac             = 1
        Float   mut_frac             = 0.5
        Float   sv_frac              = 1.0
        Float   sim_err              = 0.0
        Boolean ignore_snps          = false
        Boolean ignore_ref           = false
        Boolean ignore_sanity_check  = false
        Boolean single_ended         = false
        Boolean tag_reads            = false
        Boolean ignore_pileup        = false
        Boolean det_base_changes     = false
        Boolean keep_secondary_reads = false
        String? seed
        String  mutation_type
        String  output_basename
        String  docker_image
        Int     addldisk             = 200
        Int     mem_size             = 8
        Int     preemptible          = 1
    }

    Int bam_size        = ceil((size(input_bam, "GB") * 2.5) + size(input_bai, "GB"))
    Int fasta_size      = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB"))
    Int ref_size        = ceil(size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB"))
    Int final_disk_size = bam_size + fasta_size + ref_size + addldisk

    command <<<
        set -euxo pipefail

        mkdir --parents INPUT
        mkdir --parents OUTPUT
        input_bam_name=$(basename "~{input_bam}")
        input_bai_name=$(basename "~{input_bai}")
        mv "~{input_bam}" "INPUT/${input_bam_name}"
        mv "~{input_bai}" "INPUT/${input_bai_name}"
    
        $MGBPMBIOFXPATH/bamsurgeon/bin/create_bed.py   \
            --json-input "~{write_json(mutation_bed)}" \
            --mutation-type "~{mutation_type}"         \
            --out-bed-name "INPUT/target_regions.bed"

        $MGBPMBIOFXPATH/bamsurgeon/bin/create_sh.py         \
            --snvfrac "~{snv_frac}"                         \
            --mutfrac "~{mut_frac}"                         \
            --numsnvs "~{nums_snvs}"                        \
            --cnvfile "~{cnv_file}"                         \
            --coverdiff "~{cover_diff}"                     \
            --haplosize "~{haplo_size}"                     \
            --mindepth "~{min_depth}"                       \
            --maxdepth "~{max_depth}"                       \
            --minmutreads "~{min_mut_reads}"                \
            --avoidreads "~{avoid_reads}"                   \
            --ignoresnps "~{ignore_snps}"                   \
            --ignoreref "~{ignore_ref}"                     \
            --insane "~{ignore_sanity_check}"               \
            --single-ended "~{single_ended}"                \
            --maxopen "~{max_open_files}"                   \
            --tagreads "~{tag_reads}"                       \
            --ignorepileup "~{ignore_pileup}"               \
            --seed "~{seed}"                                \
            --det-base-changes "~{det_base_changes}"        \
            --maxlibsize "~{max_lib_size}"                  \
            --kmer-size "~{kmer_size}"                      \
            --svfrac "~{sv_frac}"                           \
            --minctglen "~{min_contig_gen}"                 \
            --donorbam "~{donor_bam}"                       \
            --ismean "~{mean_insert_size}"                  \
            --issd "~{insert_size_stdev}"                   \
            --simerr "~{sim_err}"                           \
            --inslib "~{insert_library}"                    \
            --keepsecondary "~{keep_secondary_reads}"       \
            --picardjar "$MGBPMBIOFXPATH/picard/picard.jar" \
            --aligner "mem"                                 \
            --mutation-type "~{mutation_type}"              \
            --outbam "OUTPUT/~{output_basename}"            \
            --out-script-name "INPUT/run_bamsurgeon.sh"     \
            --script-to-run "$MGBPMBIOFXPATH/bamsurgeon/bin/add~{mutation_type}.py"

        CNV_FILE=""
        [ ! -z "~{cnv_file}" ] && CNV_FILE="-c ~{cnv_file}"

        AVOID_READS=""
        [ ! -z "~{avoid_reads}" ] && AVOID_READS="-a ~{avoid_reads}"

        DONOR_BAM=""
        [ ! -z "~{donor_bam}" ] && DONOR_BAM="-d ~{donor_bam}"

        bash INPUT/run_bamsurgeon.sh \
            -v "INPUT/target_regions.bed" \
            -b "INPUT/${input_bam_name}" \
            -r "~{ref_fasta}" \
            $CNV_FILE $AVOID_READS $DONOR_BAM

        out_vcf_base=$(basename "~{output_basename}" ".bam")
        mv --target-directory=OUTPUT "${out_vcf_base}.add~{mutation_type}.target_regions.vcf"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File target_bed        = "INPUT/target_regions.bed"
        File bamsurgeon_script = "INPUT/run_bamsurgeon.sh"
        File mut_bam           = "OUTPUT/~{output_basename}"
        File mut_vcf           = "OUTPUT/" + sub(basename(output_basename), ".bam", "") + ".add~{mutation_type}.target_regions.vcf"
    }
}

task IndexBAMTask {
    input {
        File   input_bam
        String output_basename = sub(basename(input_bam), ".bam", "")
        String docker_image
        Int    addldisk        = 10
        Int    mem_size        = 4
        Int    preemptible     = 1
    }

    Int bam_size        = ceil(size(input_bam, "GB") * 3)
    Int final_disk_size = bam_size + addldisk

    command <<<
        set -euxo pipefail

        mkdir --parents OUTPUT
        samtools sort -o "OUTPUT/~{output_basename}.bam" "~{input_bam}"
        samtools index -b "OUTPUT/~{output_basename}.bam" "OUTPUT/~{output_basename}.bam.bai"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_bam = "OUTPUT/~{output_basename}.bam"
        File output_bai = "OUTPUT/~{output_basename}.bam.bai"
    }
}