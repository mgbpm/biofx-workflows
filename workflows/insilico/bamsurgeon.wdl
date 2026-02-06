version 1.0

import "../../steps/Utilities.wdl"
import "../../steps/FileUtils.wdl"
import "../../steps/IgvReport.wdl"
import "../pgxrisk/PGxWorkflow.wdl"
import "../pgxrisk/RiskAllelesWorkflow.wdl"
import "InsilicoTasks.wdl"

workflow BamsurgeonWorkflow {
    input {
        Array[MutationBED] mutation_bed
        File    basefile_idx
        File    basefile
        String  mutation_type
        String  output_basename
        File    ref_fasta
        File    ref_amb
        File    ref_ann
        File    ref_bwt
        File    ref_fai
        File    ref_pac
        File    ref_sa
        File    ref_dict
        File    dbsnp_vcf
        File    dbsnp_vcf_index
        Boolean generate_reports        = false
        String  pgx_test_code           = "lmPGX-pnlD_L"
        File    pgx_fileset             = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_20241004.tar"
        File    pgx_roi                 = "gs://lmm-reference-data/pgx/lmPGX-pnlD_L_genotyping.bed"
        String  risk_test_code          = "lmRISK-pnlB_L"
        File    risk_fileset            = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_20230105.tar"
        File    risk_roi                = "gs://lmm-reference-data/risk/lmRISK-pnlB_L_genotyping-chr_20230628.bed"
        String  bamsurgeon_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bamsurgeon:20240229"
        String  igvreport_docker_image  = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/igvreport:20230511"
        String  orchutils_docker_image  = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20251222"
        String  pgx_docker_image        = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/pgx:20241007"
        String  risk_docker_image       = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/risk:20240129"
        String  samtools_docker_image   = "biocontainers/samtools:v1.9-4-deb_cv1"
        String  gcp_project_id          = "mgb-lmm-gcp-infrast-1651079146"
        String  workspace_name
    }

    Boolean is_bam = if basename(basefile) == basename(basefile, ".bam") + ".bam" then true else false
    Boolean is_cram = if basename(basefile) == basename(basefile, ".cram") + ".cram" then true else false

    if (is_bam) {
        if (basename(basefile_idx) != (basename(basefile_idx)) + ".bai") {
            String BamIndexError = "BAM index should have BAM file name + '.bai'"
        }

        Array[String] bam_filetypes = ["bam", "bai"]
    }
    if (is_cram) {
        if (basename(basefile_idx) != (basename(basefile_idx)) + ".crai") {
            String CramIndexError = "CRAM index should have CRAM file name + '.crai'"
        }

        Array[String] cram_filetypes = ["cram", "crai"]
    }

    if ((mutation_type != "snv") && (mutation_type != "indel") && (mutation_type != "sv")) {
        String InputParameterError = "Acceptable mutation types to introduce are 'snv', 'indel' or 'sv.'"
    }
    
    String error_message = select_first([BamIndexError, CramIndexError, InputParameterError, ""])
    if (error_message != "") {
        call Utilities.FailTask as InputFailure {
            input:
                error_message = error_message
        }
    }

    if (error_message == "") {
        String file_location = sub("" + basefile, "/" + basename(basefile), "")
        Array[String] file_types = select_first([bam_filetypes, cram_filetypes])

        if (("s3" + sub(file_location, "s3", "")) == file_location) {
            call FileUtils.FetchFilesTask as FetchFiles {
                input:
                    data_location = file_location,
                    recursive = true,
                    file_types = file_types,
                    file_match_keys = [ basename(basefile, "."  + file_types[0]) ],
                    docker_image = orchutils_docker_image,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name
            }
        }

        if (is_cram) {
            call InsilicoTasks.CramToBamTask as ConvertCram {
                input:
                    input_cram = select_first([FetchFiles.cram]),
                    input_crai = select_first([FetchFiles.crai]),
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    docker_image = samtools_docker_image
            }
        }

        call InsilicoTasks.RunBamsurgeonTask as RunBamsurgeon {
            input:
                mutation_bed = mutation_bed,
                input_bam = select_first([ConvertCram.output_bam, FetchFiles.bam, basefile]),
                input_bai = select_first([ConvertCram.output_bai, FetchFiles.bai, basefile_idx]),
                ref_fasta = ref_fasta,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_bwt = ref_bwt,
                ref_fai = ref_fai,
                ref_pac = ref_pac,
                ref_sa = ref_sa,
                mutation_type = mutation_type,
                output_basename = output_basename + ".bam",
                docker_image = bamsurgeon_docker_image
        }

        call InsilicoTasks.IndexBAMTask as IndexBAM {
            input:
                input_bam = RunBamsurgeon.mut_bam,
                docker_image = samtools_docker_image
        }

        call IgvReport.IgvReportFromGenomePanelsBedTask as MutIGVReport {
            input:
                bed_file = RunBamsurgeon.target_bed,
                sample_bam = IndexBAM.output_bam,
                sample_bai = IndexBAM.output_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fai,
                output_basename = output_basename + "_IGVreport",
                docker_image = igvreport_docker_image
        }

        if (generate_reports) {
            call PGxWorkflow.PGxWorkflow as RunPGx {
                input:
                    input_cram = IndexBAM.output_bam,
                    input_crai = IndexBAM.output_bai,
                    sample_id = output_basename,
                    accession_id = output_basename,
                    test_code = pgx_test_code,
                    reference_fasta = ref_fasta,
                    reference_fasta_fai = ref_fai,
                    reference_dict = ref_dict,
                    roi_bed = pgx_roi,
                    dbsnp = dbsnp_vcf,
                    dbsnp_vcf_index = dbsnp_vcf_index,
                    workflow_fileset = pgx_fileset,
                    mgbpmbiofx_docker_image = pgx_docker_image
            }

            call RiskAllelesWorkflow.RiskAllelesWorkflow as RunRisk {
                input:
                    input_cram = IndexBAM.output_bam,
                    input_crai = IndexBAM.output_bai,
                    sample_id = output_basename,
                    accession_id = output_basename,
                    test_code = risk_test_code,
                    reference_fasta = ref_fasta,
                    reference_fasta_fai = ref_fai,
                    reference_dict = ref_dict,
                    roi_bed = risk_roi,
                    dbsnp = dbsnp_vcf,
                    dbsnp_vcf_index = dbsnp_vcf_index,
                    workflow_fileset = risk_fileset,
                    gcp_project_id = gcp_project_id,
                    workspace_name = workspace_name,
                    mgbpmbiofx_docker_image = risk_docker_image
            }
        }
    }

    output {
        File? mutated_bam = IndexBAM.output_bam
        File? mutated_bai = IndexBAM.output_bai
        File? igv_report = MutIGVReport.igv_report
        File? pgx_summary_report = RunPGx.summary_report
        File? pgx_details_report = RunPGx.details_report
        File? pgx_genotype_xlsx = RunPGx.genotype_xlsx
        File? pgx_genotype_txt = RunPGx.genotype_txt
        File? risk_report = RunRisk.risk_report
        File? risk_genotype_xlsx = RunRisk.genotype_xlsx
        File? risk_genotype_txt = RunRisk.genotype_txt
    }
}