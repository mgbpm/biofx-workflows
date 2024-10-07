version 1.0

import "../../steps/Utilities.wdl"
import "../../steps/IgvReport.wdl"
import "../pgxrisk/PGxWorkflow.wdl"
import "../pgxrisk/RiskAllelesWorkflow.wdl"

workflow BamsurgeonWorkflow {
    input {
        ## BAMSURGEON INPUTS
        # File naming inputs
        String output_files_base
        # Mutation to run
        String mutation_type
        # Input BAM files
        File input_bai_file
        File input_bam_file
        # Input for target BED file
        Array[MutationBED] mutation_bed_input
        # Reference files
        File ref_fasta
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_fai
        File ref_pac
        File ref_sa
        File? ref_dict
        File? dbsnp_vcf
        File? dbsnp_vcf_index
        # Optional inputs for all mutation types
        File? cnv_file
        Boolean tag_reads = false
        String? seed
        # Optional SNV & Indel inputs
        Float snv_frac = 1
        Float mut_frac = 0.5
        Int nums_snvs = 0
        Float cover_diff = 0.9
        Int min_depth = 10
        Int max_depth = 2000
        Int min_mut_reads = 3
        File? avoid_reads
        Boolean ignore_sanity_check = false
        Boolean single_ended = false
        Int max_open_files = 1000
        Boolean ignore_pileup = false
        # Optional SNV-only inputs
        Int haplo_size = 0
        Boolean ignore_snps = false
        Boolean ignore_ref = false
        # Optional Indel-only inputs
        Boolean det_base_changes = false
        # Optional SV-only inputs
        Int max_lib_size = 600
        Int kmer_size = 31
        Float sv_frac = 1.0
        Int min_contig_gen = 4000
        File? donor_bam
        File? donor_bai
        Int mean_insert_size = 300
        Int insert_size_stdev = 70
        Float sim_err = 0.0
        File? insert_library
        Boolean keep_secondary_reads = false
        # Bamsurgeon docker image
        String bamsurgeon_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bamsurgeon:20240229"
        # SAMtools docker image
        String samtools_docker_image = "biocontainers/samtools:v1.9-4-deb_cv1"
        # IGV docker image
        String igvreport_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511"

        ## PGX/RISK INPUTS
        # Toggle PGx/Risk
        Boolean run_pgx = true
        Boolean run_risk = true
        # PGx inputs
        String pgx_test_code = "lmPGX-pnlD_L"
        String pgx_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20241007"
        File pgx_workflow_fileset = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_20241004.tar"
        File pgx_roi_bed = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_genotyping.bed"
        # Risk inputs
        String risk_alleles_test_code = "lmRISK-pnlB_L"
        String risk_alleles_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129"
        File risk_alleles_workflow_fileset = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar"
        File risk_alleles_roi_bed = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed"
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String? workspace_name
    }

    ## Test that mutation_type is either snv, indel, or sv
    if ((mutation_type != "snv") && (mutation_type != "indel") && (mutation_type != "sv")) {
        call Utilities.FailTask as InputParameterError {
            input:
                error_message = "Acceptable mutation types to introduce are 'snv', 'indel' or 'sv.'"
        }
    }
    ## Test that the BAM index is properly named
    if (basename(input_bai_file) != (basename(input_bam_file)) + ".bai") {
        call Utilities.FailTask as BAMIndexNameError {
            input:
                error_message = "BAM index should have BAM name + '.bai'"
        }
    }
    ## Test that the inputs for PGx/Risk are included
    if (run_pgx || run_risk) {
        if (!defined(ref_dict) || !defined(dbsnp_vcf) || !defined(dbsnp_vcf_index)) {
            call Utilities.FailTask as PGxAndRiskInputsError {
                input:
                    error_message = "ref_dict, dbsnp_vcf, and dbsnp_vcf_index must be defined to run PGx and/or Risk."
            }
        }
    }
    if (run_risk && !defined(workspace_name)) {
        call Utilities.FailTask as WorkspaceNameError {
            input:
                error_message = "Workspace name must be defined to run Risk pipeline."
        }
    }

    ## Create input files for bamsurgeon
    call RunBamsurgeonTask as RunBamsurgeon {
        input:
            mutation_bed_input = mutation_bed_input,
            input_bam_file = input_bam_file,
            input_bai_file = input_bai_file,
            snv_frac = snv_frac,
            mut_frac = mut_frac,
            nums_snvs = nums_snvs,
            cnv_file = cnv_file,
            cover_diff = cover_diff,
            haplo_size = haplo_size,
            min_depth = min_depth,
            max_depth = max_depth,
            min_mut_reads = min_mut_reads,
            avoid_reads = avoid_reads,
            ignore_snps = ignore_snps,
            ignore_ref = ignore_ref,
            ignore_sanity_check = ignore_sanity_check,
            single_ended = single_ended,
            max_open_files = max_open_files,
            tag_reads = tag_reads,
            ignore_pileup = ignore_pileup,
            seed = seed,
            det_base_changes = det_base_changes,
            max_lib_size = max_lib_size,
            kmer_size = kmer_size,
            sv_frac = sv_frac,
            min_contig_gen = min_contig_gen,
            donor_bam = donor_bam,
            donor_bai = donor_bai,
            mean_insert_size = mean_insert_size,
            insert_size_stdev = insert_size_stdev,
            sim_err = sim_err,
            insert_library = insert_library,
            keep_secondary_reads = keep_secondary_reads,
            ref_fasta = ref_fasta,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_fai = ref_fai,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            mutation_type = mutation_type,
            output_bam_name = output_files_base + ".bam",
            docker_image = bamsurgeon_docker_image
    }
    call IndexBAMTask as IndexBAM {
        input:
            input_bam = RunBamsurgeon.mut_bam,
            docker_image = samtools_docker_image
    }
    ## Generate IGV report
    call IgvReport.IgvReportFromGenomePanelsBedTask as MutIGVReport {
        input:
            bed_file = RunBamsurgeon.target_bed,
            sample_bam = IndexBAM.output_bam,
            sample_bai = IndexBAM.output_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fai,
            output_basename = output_files_base + "_IGVreport",
            docker_image = igvreport_docker_image
    }
    ## Run PGx & Risk
    if (run_pgx) {
        call PGxWorkflow.PGxWorkflow as RunPGx {
            input:
                input_cram = IndexBAM.output_bam,
                input_crai = IndexBAM.output_bai,
                sample_id = output_files_base,
                accession_id = output_files_base,
                test_code = pgx_test_code,
                reference_fasta = ref_fasta,
                reference_fasta_fai = ref_fai,
                reference_dict = select_first([ref_dict]),
                roi_bed = pgx_roi_bed,
                dbsnp = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                workflow_fileset = pgx_workflow_fileset,
                mgbpmbiofx_docker_image = pgx_docker_image
        }
    }
    if (run_risk) {
        call RiskAllelesWorkflow.RiskAllelesWorkflow as RunRisk {
            input:
                input_cram = IndexBAM.output_bam,
                input_crai = IndexBAM.output_bai,
                sample_id = output_files_base,
                accession_id = output_files_base,
                test_code = risk_alleles_test_code,
                reference_fasta = ref_fasta,
                reference_fasta_fai = ref_fai,
                reference_dict = select_first([ref_dict]),
                roi_bed = risk_alleles_roi_bed,
                dbsnp = select_first([dbsnp_vcf]),
                dbsnp_vcf_index = select_first([dbsnp_vcf_index]),
                workflow_fileset = risk_alleles_workflow_fileset,
                gcp_project_id = gcp_project_id,
                workspace_name = select_first([workspace_name]),
                mgbpmbiofx_docker_image = risk_alleles_docker_image
        }
    }

    output {
        # Bamsurgeon intermediary files
        File target_bed = RunBamsurgeon.target_bed
        File bamsurgeon_script = RunBamsurgeon.bamsurgeon_script
        # Bamsurgeon outputs
        File mutated_bam = IndexBAM.output_bam
        File mutated_bai = IndexBAM.output_bai
        File mutated_vcf = RunBamsurgeon.mut_vcf
        # IGV report
        File igv_report = MutIGVReport.igv_report
        # PGx output
        File? pgx_summary_report = RunPGx.summary_report
        File? pgx_details_report = RunPGx.details_report
        File? pgx_genotype_xlsx = RunPGx.genotype_xlsx
        File? pgx_genotype_txt = RunPGx.genotype_txt
        # Risk output
        File? risk_alleles_report = RunRisk.risk_report
        File? risk_alleles_genotype_xlsx = RunRisk.genotype_xlsx
        File? risk_alleles_genotype_txt = RunRisk.genotype_txt
    }
}

struct MutationBED {
    String chr
    Int start
    Int end
    String? mut_type
    Float? vaf
    String? base
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

task RunBamsurgeonTask {
    input {
        # Input array
        Array[MutationBED] mutation_bed_input
        # Input BAM and index
        File input_bam_file
        File input_bai_file
        # Inputs to create bash and BED files
        Float snv_frac = 1
        Float mut_frac = 0.5
        Int nums_snvs = 0
        File? cnv_file
        Float cover_diff = 0.9
        Int haplo_size = 0
        Int min_depth = 10
        Int max_depth = 2000
        Int min_mut_reads = 3
        File? avoid_reads
        Boolean ignore_snps = false
        Boolean ignore_ref = false
        Boolean ignore_sanity_check = false
        Boolean single_ended = false
        Int max_open_files = 1000
        Boolean tag_reads = false
        Boolean ignore_pileup = false
        String? seed
        Boolean det_base_changes = false
        Int max_lib_size = 600
        Int kmer_size = 31
        Float sv_frac = 1.0
        Int min_contig_gen = 4000
        File? donor_bam
        File? donor_bai
        Int mean_insert_size = 300
        Int insert_size_stdev = 70
        Float sim_err = 0.0
        File? insert_library
        Boolean keep_secondary_reads = false
        # Reference files
        File ref_fasta
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_fai
        File ref_pac
        File ref_sa
        # Output-specific inputs
        String mutation_type
        String output_bam_name
        # Docker inputs
        String docker_image
        Int disk_size = ceil((size(input_bam_file, "GB") * 2.5)) + ceil(size(input_bai_file, "GB")) + ceil(size(ref_fasta, "GB")) + ceil(size(ref_amb, "GB")) + ceil(size(ref_ann, "GB")) + ceil(size(ref_bwt, "GB")) + ceil(size(ref_fai, "GB")) + ceil(size(ref_pac, "GB")) + ceil(size(ref_sa, "GB")) + 10
        Int mem_size = 8
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Make input & output directories
        mkdir --parents INPUT
        mkdir --parents OUTPUT

        # Move BAM and index to input dir
        input_bam_name=$(basename "~{input_bam_file}")
        input_bai_name=$(basename "~{input_bai_file}")
        mv "~{input_bam_file}" "INPUT/${input_bam_name}"
        mv "~{input_bai_file}" "INPUT/${input_bai_name}"
    
        # Create the BED file for bamsurgeon to use
        $MGBPMBIOFXPATH/bamsurgeon/bin/create_bed.py \
            --json-input "~{write_json(mutation_bed_input)}" \
            --mutation-type "~{mutation_type}" \
            --out-bed-name "INPUT/target_regions.bed"

        # Create bash script that will run bamsurgeon
        $MGBPMBIOFXPATH/bamsurgeon/bin/create_sh.py \
            --snvfrac "~{snv_frac}" \
            --mutfrac "~{mut_frac}" \
            --numsnvs "~{nums_snvs}" \
            --cnvfile "~{cnv_file}" \
            --coverdiff "~{cover_diff}" \
            --haplosize "~{haplo_size}" \
            --mindepth "~{min_depth}" \
            --maxdepth "~{max_depth}" \
            --minmutreads "~{min_mut_reads}" \
            --avoidreads "~{avoid_reads}" \
            --ignoresnps "~{ignore_snps}" \
            --ignoreref "~{ignore_ref}" \
            --insane "~{ignore_sanity_check}" \
            --single-ended "~{single_ended}" \
            --maxopen "~{max_open_files}" \
            --tagreads "~{tag_reads}" \
            --ignorepileup "~{ignore_pileup}" \
            --seed "~{seed}" \
            --det-base-changes "~{det_base_changes}" \
            --maxlibsize "~{max_lib_size}" \
            --kmer-size "~{kmer_size}" \
            --svfrac "~{sv_frac}" \
            --minctglen "~{min_contig_gen}" \
            --donorbam "~{donor_bam}" \
            --ismean "~{mean_insert_size}" \
            --issd "~{insert_size_stdev}" \
            --simerr "~{sim_err}" \
            --inslib "~{insert_library}" \
            --keepsecondary "~{keep_secondary_reads}" \
            --picardjar "$MGBPMBIOFXPATH/picard/picard.jar" \
            --aligner "mem" \
            --mutation-type "~{mutation_type}" \
            --outbam "OUTPUT/~{output_bam_name}" \
            --out-script-name "INPUT/run_bamsurgeon.sh" \
            --script-to-run "$MGBPMBIOFXPATH/bamsurgeon/bin/add~{mutation_type}.py"

        # Get a list of the optional parameters being run
        CNV_FILE=""
        [ ! -z "~{cnv_file}" ] && CNV_FILE="-c ~{cnv_file}"

        AVOID_READS=""
        [ ! -z "~{avoid_reads}" ] && AVOID_READS="-a ~{avoid_reads}"

        DONOR_BAM=""
        [ ! -z "~{donor_bam}" ] && DONOR_BAM="-d ~{donor_bam}"

        OPTIONAL_PARAMS="$CNV_FILE $AVOID_READS $DONOR_BAM"

        # Run the bamsurgeon bash script with the appropriate parameters
        bash INPUT/run_bamsurgeon.sh \
            -v "INPUT/target_regions.bed" \
            -b "INPUT/${input_bam_name}" \
            -r "~{ref_fasta}" \
            $OPTIONAL_PARAMS

        # Move output VCF to OUTPUT dir
        out_vcf_base=$(basename "~{output_bam_name}" ".bam")
        mv --target-directory=OUTPUT "${out_vcf_base}.add~{mutation_type}.target_regions.vcf"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        # Intermediary files
        File target_bed = "INPUT/target_regions.bed"
        File bamsurgeon_script = "INPUT/run_bamsurgeon.sh"
        # Final outputs from bamsurgeon
        File mut_bam = "OUTPUT/~{output_bam_name}"
        File mut_vcf = "OUTPUT/" + sub(basename(output_bam_name), ".bam", "") + ".add~{mutation_type}.target_regions.vcf"
    }
}

task IndexBAMTask {
    input {
        File input_bam
        String output_basename = sub(basename(input_bam), ".bam", "") + ".sorted"
        String docker_image
        Int disk_size = ceil((size(input_bam, "GB") * 3) + 10)
        Int mem_size = 4
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        mkdir --parents OUTPUT

        # Sort the input BAM
        samtools sort -o "OUTPUT/~{output_basename}.bam" "~{input_bam}"

        # Index the output mutated BAM
        samtools index -b "OUTPUT/~{output_basename}.bam" "OUTPUT/~{output_basename}.bam.bai"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_bam = "OUTPUT/~{output_basename}.bam"
        File output_bai = "OUTPUT/~{output_basename}.bam.bai"
    }
}