version 1.0

import "../../steps/Utilities.wdl"
import "../../steps/IgvReport.wdl"
import "../pgxrisk/PGxWorkflow.wdl"
import "../pgxrisk/RiskAllelesWorkflow.wdl"

workflow insiMWorkflow {
    input {
        ## INSIM INPUTS
        # Input BAM files
        File input_bam
        File input_bai
        # Reference files
        File ref_fasta
        File ref_fai
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File? ref_dict
        File? dbsnp_vcf
        File? dbsnp_vcf_index
        # Target inputs
        Array[MutationTargets] targets
        Array[AmpliconBed] ampbed = []
        # insiM inputs
        String assay_type
        String mutation_type
        String? vaf
        String? inslen
        String? insseq
        String? amp_read_length
        # File naming inputs
        String sample_name
        String mutation_name
        # insiM docker image
        String insim_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/insim:testing"
        # Picard docker image
        String picard_docker_image = "biocontainers/picard:v1.139_cv3"
        # BWA docker image
        String bwa_docker_image = "biocontainers/bwa:v0.7.17_cv1"
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
        String pgx_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20240614"
        File pgx_workflow_fileset = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_20240606.tar"
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

    ## Run tests on inputs to ensure their compatibility
    # Test for assay type
    if ((assay_type != "capture") && (assay_type != "amplicon")) {
        String AssayTypeFail = "Assay type should either be 'capture' or 'amplicon'."
    }
    # Test for mutation type
    if ((mutation_type != "snv") && (mutation_type != "ins") && (mutation_type != "del") && (mutation_type != "dup") && (mutation_type != "mix")) {
        String MutationTypeFail = "Acceptable mutation types are 'snv', 'ins', 'del', 'dup', or 'mix.'"
    }
    # Test for VAF if mutation type is not "mix"
    if ((mutation_type != "mix") && (!defined(vaf))) {
        String VAFFail = "VAF input is required for non-'mix' (snv, ins, del, dup) mutations."
    }
    # Test for input ampbed_array and amp_read_length when input is an amplicon assay
    if ((assay_type == "amplicon") && (!defined(amp_read_length))) {
        String AmpReadLengthFail = "Amplicon read length not defined"
    }
    # Fail the workflow if there is an error
    String input_error_message = select_first([AssayTypeFail, MutationTypeFail, VAFFail, AmpReadLengthFail, ""])
    if (input_error_message != "") {
        call Utilities.FailTask as InputParameterError {
            input:
                error_message = input_error_message
        }
    }

    ## Mutate input BAM based on assay type
    String output_basename = sample_name + "_" + mutation_name + ".bam"
    if (assay_type == 'amplicon') {
        call MutateAmpliconAssayTask as MutateAmplicon {
            input:
                input_bam = input_bam,
                input_bai = input_bai,
                assay_type = assay_type,
                targets = targets,
                mutation_type = mutation_type,
                ref_fasta = ref_fasta,
                amp_read_length = amp_read_length,
                vaf = vaf,
                inslen = inslen,
                insseq = insseq,
                ampbed = ampbed,
                output_basename = output_basename,
                docker_image = insim_docker_image
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
                ref_fasta = ref_fasta,
                vaf = vaf,
                inslen = inslen,
                insseq = insseq,
                output_basename = output_basename,
                docker_image = insim_docker_image
        }
    }

    ## Convert mutations to BAM
    # Convert FASTQ to SAM
    call ConvertToSAMTask as ConvertToSAM {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            in_fastq1 = select_first([MutateCapture.output_fastq1, MutateAmplicon.output_fastq1]),
            in_fastq2 = select_first([MutateCapture.output_fastq2, MutateAmplicon.output_fastq2]),
            output_basename = output_basename,
            docker_image = bwa_docker_image
    }
    # Convert SAM to BAM
    call ConvertToBAMTask as ConvertToBAM {
        input:
            input_sam = ConvertToSAM.output_sam,
            output_basename = output_basename,
            docker_image = picard_docker_image
    }
    # Index final BAM
    call IndexBAMTask as IndexBAM {
        input:
            output_basename = output_basename,
            input_bam = ConvertToBAM.output_dedup_bam,
            docker_image = samtools_docker_image
    }

    ## Generate IGV report
    call IgvReport.IgvReportFromGenomePanelsBedTask as MutIGVReport {
        input:
            bed_file = select_first([MutateAmplicon.target_bed, MutateCapture.target_bed]),
            sample_bam = ConvertToBAM.output_dedup_bam,
            sample_bai = IndexBAM.output_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fai,
            output_basename = sample_name + "_" + mutation_name + "_IGVreport",
            docker_image = igvreport_docker_image
    }

     ## Run PGx & Risk
    if (run_pgx) {
        call PGxWorkflow.PGxWorkflow as MutBamPGx {
            input:
                input_cram = ConvertToBAM.output_dedup_bam,
                input_crai = IndexBAM.output_bai,
                sample_id = sample_name + "_" + mutation_name,
                accession_id = sample_name + "_" + mutation_name,
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
        call RiskAllelesWorkflow.RiskAllelesWorkflow as MutBamRisk {
            input:
                input_cram = ConvertToBAM.output_dedup_bam,
                input_crai = IndexBAM.output_bai,
                sample_id = sample_name + "_" + mutation_name,
                accession_id = sample_name + "_" + mutation_name,
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
        # Files from converting to BAM
        File output_dedup_bam = ConvertToBAM.output_dedup_bam
        File output_dedup_metrics = ConvertToBAM.metrics_file
        File output_dedup_bai = IndexBAM.output_bai
        # IGV report
        File mut_igv_report = MutIGVReport.igv_report
        # PGx output
        File? pgx_summary_report = MutBamPGx.summary_report
        File? pgx_details_report = MutBamPGx.details_report
        File? pgx_genotype_xlsx = MutBamPGx.genotype_xlsx
        File? pgx_genotype_txt = MutBamPGx.genotype_txt
        # Risk output
        File? risk_alleles_report = MutBamRisk.risk_report
        File? risk_alleles_genotype_xlsx = MutBamRisk.genotype_xlsx
        File? risk_alleles_genotype_txt = MutBamRisk.genotype_txt
    }
}

struct MutationTargets {
    String chrom
    Int start
    Int end
    String type
    Float vaf
    String? ins_seq
    Int? ins_len
}

struct AmpliconBed {
    String chrom
    Int start
    Int end
}

task MutateAmpliconAssayTask {
    input {
        # Input BAM
        File input_bam
        File input_bai
        # insiM inputs
        String assay_type
        Array[MutationTargets] targets
        Array[AmpliconBed] ampbed = []
        String mutation_type
        File ref_fasta
        String vaf = ""
        String inslen = ""
        String insseq = ""
        String? amp_read_length
        # Docker inputs
        String output_basename
        String docker_image
        Int disk_size = ceil((size(input_bam, "GB") * 3)) + ceil(size(ref_fasta, "GB")) + 10
        Int mem_size = 16
    }

    command <<<
        set -euxo pipefail

        # Create config file
        $MGBPMBIOFXPATH/insiM/create_insiM_input.py \
            --assay "~{assay_type}" \
            --input-bam "~{input_bam}" \
            --target-json-file "~{write_json(targets)}" \
            --mutation-type "~{mutation_type}" \
            --genome-ref "~{ref_fasta}" \
            --config-file "amplicon_config.txt" \
            --bed-name "amplicon_targets.bed" \
            --ampliconbed-json-file "~{write_json(ampbed)}" \
            --read-length "~{amp_read_length}" \
            --vaf "~{vaf}" \
            --inslen "~{inslen}" \
            --insseq "~{insseq}" \
            --output-file "~{output_basename}"

        # Run insiM with the config file
        $MGBPMBIOFXPATH/insiM/insiM_v2.0.py -config "amplicon_config.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
    }

    output {
        # Intermediary files
        File insiM_config = "amplicon_config.txt"
        File target_bed = "amplicon_targets.bed"
        # Final insiM outputs
        File output_fastq1 = "~{output_basename}_R1_001.fastq"
        File output_fastq2 = "~{output_basename}_R2_001.fastq"
        File output_vcf = "~{output_basename}.insiM.vcf"
    }
}

task MutateCaptureAssayTask {
    input {
        # Input BAM
        File input_bam
        File input_bai
        # insiM inputs
        String assay_type
        Array[MutationTargets] targets
        String mutation_type
        File ref_fasta
        String vaf = ""
        String inslen = ""
        String insseq = ""
        # Docker inputs
        String output_basename
        String docker_image
        Int disk_size = ceil((size(input_bam, "GB") * 3)) + ceil(size(ref_fasta, "GB")) + 10
        Int mem_size = 16
    }

    command <<<
        set -euxo pipefail

        # Create config file
        $MGBPMBIOFXPATH/insiM/bin/create_insiM_input.py \
            --assay "~{assay_type}" \
            --input-bam "~{input_bam}" \
            --target-json-file "~{write_json(targets)}" \
            --mutation-type "~{mutation_type}" \
            --genome-ref "~{ref_fasta}" \
            --config-file "capture_config.txt" \
            --bed-name "capture_targets.bed" \
            --vaf "~{vaf}" \
            --inslen "~{inslen}" \
            --insseq "~{insseq}" \
            --output-file "~{output_basename}"

        # Run insiM with the config file
        $MGBPMBIOFXPATH/insiM/bin/insiM_v2.0.py -config "capture_config.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
    }

    output {
        # Intermediary files
        File insiM_config = "capture_config.txt"
        File target_bed = "capture_targets.bed"
        # Final insiM outputs
        File output_fastq1 = "~{output_basename}_R1_001.fastq"
        File output_fastq2 = "~{output_basename}_R2_001.fastq"
        File output_vcf = "~{output_basename}.insiM.vcf"
    }
}

task ConvertToSAMTask {
    input {
        # Input FASTQs
        File in_fastq1
        File in_fastq2
        # Reference files
        File ref_fasta
        File ref_fai
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        # Docker inputs
        String output_basename
        String docker_image
        Int disk_size = ceil(size(ref_fasta, "GB")) + ceil(size(ref_amb, "GB")) + ceil(size(ref_ann, "GB")) + ceil(size(ref_bwt, "GB")) + ceil(size(ref_fai, "GB")) + ceil(size(ref_pac, "GB")) + ceil(size(ref_sa, "GB")) + ceil((size(in_fastq1, "GB") * 4)) + 10
        Int mem_size = 8
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Align the FASTQ's to reference and produce SAM output
        bwa mem -M "~{ref_fasta}" "~{in_fastq1}" "~{in_fastq2}" > "~{output_basename}_aligned_reads.sam"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_sam =  "~{output_basename}_aligned_reads.sam"
    }
}

task ConvertToBAMTask{
    input {
        File input_sam
        String output_basename
        String docker_image
        String picard_path = "/opt/conda/bin/"
        Int disk_size = ceil((size(input_sam, "GB") * 2) + 10)
        Int mem_size = 4
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Sort SAM and conver to BAM
        java -jar "~{picard_path}/picard.jar" SortSam \
            INPUT="~{input_sam}" \
            OUTPUT="~{output_basename}.sorted.bam" \
            SORT_ORDER=coordinate

        # Mark and remove duplicates in BAM
        java -jar "~{picard_path}/picard.jar" MarkDuplicates \
            INPUT="~{output_basename}.sorted.bam" \
            OUTPUT="~{output_basename}.dedup.bam" \
            REMOVE_DUPLICATES=true \
            METRICS_FILE="~{output_basename}.metrics.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_sorted_bam = "~{output_basename}.sorted.bam"
        File output_dedup_bam = "~{output_basename}.dedup.bam"
        File metrics_file = "~{output_basename}.metrics.txt"
    }
}

task IndexBAMTask {
    input {
        File input_bam
        String output_basename = basename(input_bam)
        String docker_image
        Int disk_size = ceil((size(input_bam, "GB") * 1.5) + 10)
        Int mem_size = 4
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Index the output mutated BAM
        samtools index -b "~{input_bam}" "~{output_basename}.bai"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_bai = "~{output_basename}.bai"
    }
}