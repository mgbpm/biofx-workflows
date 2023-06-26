version 1.0

workflow BsurgeonWorkflow {
    input {
        # mutation to run is either snv, indel, or sv
        String mutation_to_run
        String out_bam_name
        String? seed
        String mgbpmbiofx_docker_image
        Array[SNVInput] snv_input = []
        Array[IndelInput] indel_input = []
        Array[SVInput] sv_input = []
        File bam_file
        File bai_file
        File ref_fasta
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_fai
        File ref_fasta_pac
        File ref_fasta_sa
        File? cnv_file
        File? avoid_reads
        File? donor_bam
        File? donor_bai
        File? insert_library
        Float snv_frac = 1
        Float mut_frac = 0.5
        Float sv_frac = 1.0
        Float cover_diff = 0.9
        Float sim_err = 0.0
        Int nums_snvs = 0
        Int haplo_size = 0
        Int min_depth = 10
        Int max_depth = 2000
        Int min_mut_reads = 3
        Int max_open_files = 1000
        Int max_lib_size = 600
        Int kmer_size = 31
        Int min_contig_gen = 4000
        Int mean_insert_size = 300
        Int insert_size_stdev = 70
        Boolean ignore_snps = false
        Boolean ignore_ref = false
        Boolean ignore_sanity_check = false
        Boolean single_ended = false
        Boolean tag_reads = false
        Boolean ignore_pileup = false
        Boolean deterministic_base_changes = false
        Boolean keep_secondary_reads = false
    }

    # run bamsurgeon either for snv, indel, or sv mutation
    if (mutation_to_run == "snv") {
        call CreateSNVFilesTask as CreateSNVFiles {
            input:
                mutation_to_run = mutation_to_run,
                out_bam_name = out_bam_name,
                seed = seed,
                snv_input = snv_input,
                cnv_file = cnv_file,
                avoid_reads = avoid_reads,
                donor_bam = donor_bam,
                insert_library = insert_library,
                snv_frac = snv_frac,
                mut_frac = mut_frac,
                sv_frac = sv_frac,
                cover_diff = cover_diff,
                sim_err = sim_err,
                nums_snvs = nums_snvs,
                haplo_size = haplo_size,
                min_depth = min_depth,
                max_depth = max_depth,
                min_mut_reads = min_mut_reads,
                max_open_files = max_open_files,
                max_lib_size = max_lib_size,
                kmer_size = kmer_size,
                min_contig_gen = min_contig_gen,
                mean_insert_size = mean_insert_size,
                insert_size_stdev = insert_size_stdev,
                ignore_snps = ignore_snps,
                ignore_ref = ignore_ref,
                ignore_sanity_check = ignore_sanity_check,
                single_ended = single_ended,
                tag_reads = tag_reads,
                ignore_pileup = ignore_pileup,
                deterministic_base_changes = deterministic_base_changes,
                keep_secondary_reads = keep_secondary_reads,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }

        call RunBamsurgeonTask as RunSNV {
            input:
                out_bam_name = out_bam_name,
                bam_file = bam_file,
                bai_file = bai_file,
                ref_fasta = ref_fasta,
                ref_fasta_amb = ref_fasta_amb,
                ref_fasta_ann = ref_fasta_ann,
                ref_fasta_bwt = ref_fasta_bwt,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_pac = ref_fasta_pac,
                ref_fasta_sa = ref_fasta_sa,
                target_bed = CreateSNVFiles.target_bed,
                bsurgeon_sh = CreateSNVFiles.bsurgeon_sh,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }
    }

    if (mutation_to_run == "indel"){
        call CreateIndelFilesTask as CreateIndelFiles {
            input:
                mutation_to_run = mutation_to_run,
                out_bam_name = out_bam_name,
                seed = seed,
                indel_input = indel_input,
                cnv_file = cnv_file,
                avoid_reads = avoid_reads,
                donor_bam = donor_bam,
                insert_library = insert_library,
                snv_frac = snv_frac,
                mut_frac = mut_frac,
                sv_frac = sv_frac,
                cover_diff = cover_diff,
                sim_err = sim_err,
                nums_snvs = nums_snvs,
                haplo_size = haplo_size,
                min_depth = min_depth,
                max_depth = max_depth,
                min_mut_reads = min_mut_reads,
                max_open_files = max_open_files,
                max_lib_size = max_lib_size,
                kmer_size = kmer_size,
                min_contig_gen = min_contig_gen,
                mean_insert_size = mean_insert_size,
                insert_size_stdev = insert_size_stdev,
                ignore_snps = ignore_snps,
                ignore_ref = ignore_ref,
                ignore_sanity_check = ignore_sanity_check,
                single_ended = single_ended,
                tag_reads = tag_reads,
                ignore_pileup = ignore_pileup,
                deterministic_base_changes = deterministic_base_changes,
                keep_secondary_reads = keep_secondary_reads,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }

        call RunBamsurgeonTask as RunIndel {
            input:
                out_bam_name = out_bam_name,
                bam_file = bam_file,
                bai_file = bai_file,
                ref_fasta = ref_fasta,
                ref_fasta_amb = ref_fasta_amb,
                ref_fasta_ann = ref_fasta_ann,
                ref_fasta_bwt = ref_fasta_bwt,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_pac = ref_fasta_pac,
                ref_fasta_sa = ref_fasta_sa,
                target_bed = CreateIndelFiles.target_bed,
                bsurgeon_sh = CreateIndelFiles.bsurgeon_sh,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }
    }

    if (mutation_to_run == "sv"){
        call CreateSVFilesTask as CreateSVFiles {
            input:
                mutation_to_run = mutation_to_run,
                out_bam_name = out_bam_name,
                seed = seed,
                sv_input = sv_input,
                cnv_file = cnv_file,
                avoid_reads = avoid_reads,
                donor_bam = donor_bam,
                donor_bai = donor_bai,
                insert_library = insert_library,
                snv_frac = snv_frac,
                mut_frac = mut_frac,
                sv_frac = sv_frac,
                cover_diff = cover_diff,
                sim_err = sim_err,
                nums_snvs = nums_snvs,
                haplo_size = haplo_size,
                min_depth = min_depth,
                max_depth = max_depth,
                min_mut_reads = min_mut_reads,
                max_open_files = max_open_files,
                max_lib_size = max_lib_size,
                kmer_size = kmer_size,
                min_contig_gen = min_contig_gen,
                mean_insert_size = mean_insert_size,
                insert_size_stdev = insert_size_stdev,
                ignore_snps = ignore_snps,
                ignore_ref = ignore_ref,
                ignore_sanity_check = ignore_sanity_check,
                single_ended = single_ended,
                tag_reads = tag_reads,
                ignore_pileup = ignore_pileup,
                deterministic_base_changes = deterministic_base_changes,
                keep_secondary_reads = keep_secondary_reads,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }

        call RunBamsurgeonTask as RunSV {
            input:
                out_bam_name = out_bam_name,
                bam_file = bam_file,
                bai_file = bai_file,
                ref_fasta = ref_fasta,
                ref_fasta_amb = ref_fasta_amb,
                ref_fasta_ann = ref_fasta_ann,
                ref_fasta_bwt = ref_fasta_bwt,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_pac = ref_fasta_pac,
                ref_fasta_sa = ref_fasta_sa,
                target_bed = CreateSVFiles.target_bed,
                bsurgeon_sh = CreateSVFiles.bsurgeon_sh,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }
    }

    output {
        File final_bam = select_first([RunSNV.mut_bam, RunIndel.mut_bam, RunSV.mut_bam])
    }
}

# array structures to create bed file by mutation type
struct SNVInput {
    String chr
    Int start
    Int end
    Float? vaf
    String? base
}

struct IndelInput {
    String chr
    Int start
    Int end
    Float vaf
    String mut_type
    String? insert_seq
}

struct SVInput {
    String chr
    Int start
    Int end
    String mut_type
    Float vaf
    String? insert_seq
    String? target_site_len
    String? cut_site_motif
    String? trans_chr
    Int? trans_start
    String? trans_end
    String? trans_on_chr
    String? trans_from_chr
    Int? num_dupes
}

task CreateSNVFilesTask {
    input {
        String mutation_to_run
        String out_bam_name
        String mgbpmbiofx_docker_image
        String? seed
        Array[SNVInput] snv_input = []
        File? cnv_file
        File? avoid_reads
        File? donor_bam
        File? insert_library
        Float snv_frac = 1
        Float mut_frac = 0.5
        Float sv_frac = 1.0
        Float cover_diff = 0.9
        Float sim_err = 0.0
        Int nums_snvs = 0
        Int haplo_size = 0
        Int min_depth = 10
        Int max_depth = 2000
        Int min_mut_reads = 3
        Int max_open_files = 1000
        Int max_lib_size = 600
        Int kmer_size = 31
        Int min_contig_gen = 4000
        Int mean_insert_size = 300
        Int insert_size_stdev = 70
        Boolean ignore_snps = false
        Boolean ignore_ref = false
        Boolean ignore_sanity_check = false
        Boolean single_ended = false
        Boolean tag_reads = false
        Boolean ignore_pileup = false
        Boolean deterministic_base_changes = false
        Boolean keep_secondary_reads = false
    }

    command <<<
        set -euxo pipefail

        # create target regions bed file for bamsurgeon input
            # bed file to create == "target_regions.bed"
        $MGBPMBIOFXPATH/bamsurgeon/create_bed.py \
            --json-input "~{write_json(snv_input)}" \
            --mutation-type "~{mutation_to_run}" \
            --out-bed-name "target_regions.bed"

        # create bash script that will run bamsurgeon
            # bash script to create == "run_bamsurgeon.sh"
            # bamsurgeon script to run == "addsnv.py"
        $MGBPMBIOFXPATH/bamsurgeon/create_sh.py \
            --mutation-type "~{mutation_to_run}" \
            --outbam "~{out_bam_name}" \
            --snvfrac "~{snv_frac}" \
            --mutfrac "~{mut_frac}" \
            --numsnvs "~{nums_snvs}" \
            --cnvfile "~{cnv_file}" \
            --coverdiff "~{cover_diff}" \
            --haplosize "~{haplo_size}" \
            --picardjar "/picard.jar" \
            --mindepth "~{min_depth}" \
            --maxdepth "~{max_depth}" \
            --minmutreads "~{min_mut_reads}" \
            --avoidreads "~{avoid_reads}" \
            --det "~{deterministic_base_changes}" \
            --ignoresnps "~{ignore_snps}" \
            --ignoreref "~{ignore_ref}" \
            --insane "~{ignore_sanity_check}" \
            --single "~{single_ended}" \
            --maxopen "~{max_open_files}" \
            --tagreads "~{tag_reads}" \
            --ignorepileup "~{ignore_pileup}" \
            --aligner "mem" \
            --seed "~{seed}" \
            --maxlibsize "~{max_lib_size}" \
            --kmer "~{kmer_size}" \
            --svfrac "~{sv_frac}" \
            --minctglen "~{min_contig_gen}" \
            --donorbam "~{donor_bam}" \
            --ismean "~{mean_insert_size}" \
            --issd "~{insert_size_stdev}" \
            --simerr "~{sim_err}" \
            --inslib "~{insert_library}" \
            --keepsecondary "~{keep_secondary_reads}"
            --out-script-name "run_bamsurgeon.sh"
            --script-to-run "$MGBPMBIOFXPATH/bamsurgeon/bin/addsnv.py"
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
    }

    output {
        File target_bed = "target_regions.bed"
        File bsurgeon_sh = "run_bamsurgeon.sh"
    }
}

task CreateIndelFilesTask {
    input {
        String mutation_to_run
        String out_bam_name
        String mgbpmbiofx_docker_image
        String? seed
        Array[IndelInput] indel_input = []
        File? cnv_file
        File? avoid_reads
        File? donor_bam
        File? insert_library
        Float snv_frac = 1
        Float mut_frac = 0.5
        Float sv_frac = 1.0
        Float cover_diff = 0.9
        Float sim_err = 0.0
        Int nums_snvs = 0
        Int haplo_size = 0
        Int min_depth = 10
        Int max_depth = 2000
        Int min_mut_reads = 3
        Int max_open_files = 1000
        Int max_lib_size = 600
        Int kmer_size = 31
        Int min_contig_gen = 4000
        Int mean_insert_size = 300
        Int insert_size_stdev = 70
        Boolean ignore_snps = false
        Boolean ignore_ref = false
        Boolean ignore_sanity_check = false
        Boolean single_ended = false
        Boolean tag_reads = false
        Boolean ignore_pileup = false
        Boolean deterministic_base_changes = false
        Boolean keep_secondary_reads = false
    }

    command <<<
        set -euxo pipefail

        # create target regions bed file for bamsurgeon input
            # bed file to create == "target_regions.bed"
        $MGBPMBIOFXPATH/bamsurgeon/create_bed.py \
            --json-input "~{write_json(indel_input)}" \
            --mutation-type "~{mutation_to_run}" \
            --out-bed-name "target_regions.bed"

        # create bash script that will run bamsurgeon
            # bash script to create == "run_bamsurgeon.sh"
            # bamsurgeon script to run == "addindel.py"
        $MGBPMBIOFXPATH/bamsurgeon/create_sh.py \
            --mutation-type "~{mutation_to_run}" \
            --outbam "~{out_bam_name}" \
            --snvfrac "~{snv_frac}" \
            --mutfrac "~{mut_frac}" \
            --numsnvs "~{nums_snvs}" \
            --cnvfile "~{cnv_file}" \
            --coverdiff "~{cover_diff}" \
            --haplosize "~{haplo_size}" \
            --picardjar "/picard.jar" \
            --mindepth "~{min_depth}" \
            --maxdepth "~{max_depth}" \
            --minmutreads "~{min_mut_reads}" \
            --avoidreads "~{avoid_reads}" \
            --det "~{deterministic_base_changes}" \
            --ignoresnps "~{ignore_snps}" \
            --ignoreref "~{ignore_ref}" \
            --insane "~{ignore_sanity_check}" \
            --single "~{single_ended}" \
            --maxopen "~{max_open_files}" \
            --tagreads "~{tag_reads}" \
            --ignorepileup "~{ignore_pileup}" \
            --aligner "mem" \
            --seed "~{seed}" \
            --maxlibsize "~{max_lib_size}" \
            --kmer "~{kmer_size}" \
            --svfrac "~{sv_frac}" \
            --minctglen "~{min_contig_gen}" \
            --donorbam "~{donor_bam}" \
            --ismean "~{mean_insert_size}" \
            --issd "~{insert_size_stdev}" \
            --simerr "~{sim_err}" \
            --inslib "~{insert_library}" \
            --keepsecondary "~{keep_secondary_reads}"
            --out-script-name "run_bamsurgeon.sh"
            --script-to-run "$MGBPMBIOFXPATH/bamsurgeon/bin/addindel.py"
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
    }

    output {
        File target_bed = "target_regions.bed"
        File bsurgeon_sh = "run_bamsurgeon.sh"
    }
}

task CreateSVFilesTask {
    input {
        String mutation_to_run
        String out_bam_name
        String mgbpmbiofx_docker_image
        String? seed
        Array[SVInput] sv_input = []
        File? cnv_file
        File? avoid_reads
        File? donor_bam
        File? donor_bai
        File? insert_library
        Float snv_frac = 1
        Float mut_frac = 0.5
        Float sv_frac = 1.0
        Float cover_diff = 0.9
        Float sim_err = 0.0
        Int nums_snvs = 0
        Int haplo_size = 0
        Int min_depth = 10
        Int max_depth = 2000
        Int min_mut_reads = 3
        Int max_open_files = 1000
        Int max_lib_size = 600
        Int kmer_size = 31
        Int min_contig_gen = 4000
        Int mean_insert_size = 300
        Int insert_size_stdev = 70
        Boolean ignore_snps = false
        Boolean ignore_ref = false
        Boolean ignore_sanity_check = false
        Boolean single_ended = false
        Boolean tag_reads = false
        Boolean ignore_pileup = false
        Boolean deterministic_base_changes = false
        Boolean keep_secondary_reads = false
    }

    command <<<
        set -euxo pipefail

        # create target regions bed file for bamsurgeon input
            # bed file to create == "target_regions.bed"
        $MGBPMBIOFXPATH/bamsurgeon/create_bed.py \
            --json-input "~{write_json(sv_input)}" \
            --mutation-type "~{mutation_to_run}" \
            --out-bed-name "target_regions.bed"

        # create bash script that will run bamsurgeon
            # bash script to create == "run_bamsurgeon.sh"
            # bamsurgeon script to run == "addsv.py"
        $MGBPMBIOFXPATH/bamsurgeon/create_sh.py \
            --mutation-type "~{mutation_to_run}" \
            --outbam "~{out_bam_name}" \
            --snvfrac "~{snv_frac}" \
            --mutfrac "~{mut_frac}" \
            --numsnvs "~{nums_snvs}" \
            --cnvfile "~{cnv_file}" \
            --coverdiff "~{cover_diff}" \
            --haplosize "~{haplo_size}" \
            --picardjar "/picard.jar" \
            --mindepth "~{min_depth}" \
            --maxdepth "~{max_depth}" \
            --minmutreads "~{min_mut_reads}" \
            --avoidreads "~{avoid_reads}" \
            --det "~{deterministic_base_changes}" \
            --ignoresnps "~{ignore_snps}" \
            --ignoreref "~{ignore_ref}" \
            --insane "~{ignore_sanity_check}" \
            --single "~{single_ended}" \
            --maxopen "~{max_open_files}" \
            --tagreads "~{tag_reads}" \
            --ignorepileup "~{ignore_pileup}" \
            --aligner "mem" \
            --seed "~{seed}" \
            --maxlibsize "~{max_lib_size}" \
            --kmer "~{kmer_size}" \
            --svfrac "~{sv_frac}" \
            --minctglen "~{min_contig_gen}" \
            --donorbam "~{donor_bam}" \
            --ismean "~{mean_insert_size}" \
            --issd "~{insert_size_stdev}" \
            --simerr "~{sim_err}" \
            --inslib "~{insert_library}" \
            --keepsecondary "~{keep_secondary_reads}"
            --out-script-name "run_bamsurgeon.sh"
            --script-to-run "$MGBPMBIOFXPATH/bamsurgeon/bin/addsv.py"
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
    }

    output {
        File target_bed = "target_regions.bed"
        File bsurgeon_sh = "run_bamsurgeon.sh"
    }
}

task RunBamsurgeonTask {
    input {
        File target_bed
        File bsurgeon_sh
        File bam_file
        File bai_file
        File ref_fasta
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_fai
        File ref_fasta_pac
        File ref_fasta_sa
        File? cnv_file
        File? avoid_reads
        File? donor_bam
        File? donor_bai
        File? insert_library
        String out_bam_name
        String mgbpmbiofx_docker_image
    }

    command <<<
        set -euxo pipefail

        # get a list of the optional parameters being run
        CNV_FILE=""
        [ ! -z "~{cnv_file}" ] && CNV_FILE="-c ~{cnv_file}"

        AVOID_READS=""
        [ ! -z "~{avoid_reads}" ] && AVOID_READS="-a ~{avoid_reads}"

        DONOR_BAM=""
        [ ! -z "~{donor_bam}" ] && DONOR_BAM="-d ~{donor_bam}"

        OPTIONAL_PARAMS="$CNV_FILE $AVOID_READS $DONOR_BAM"

        # run the bamsurgeon bash script with the appropriate parameters
        "~{bsurgeon_sh}" \
            -v "~{target_bed}" \
            -b "~{bam_file}" \
            -r "~{ref_fasta}" \
            $OPTIONAL_PARAMS
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "local-disk " + ((size(bam_file, "GB") * 2.5) + 20) + " HDD"
    }

    output {
        File mut_bam = out_bam_name
    }
}