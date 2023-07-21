version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"

workflow BgwgsPanelWorkflow {

    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id
        String workspace_name
        # Orchestration utils docker
        String orchutils_docker_image
        # subject, sample id and data location
        String subject_id
        String sample_id
        String sample_run_id
        String sample_lmm_id
        String sample_barcode
        String batch_id
        String sample_data_location
        String test_code
        Int fetch_disk_size = 75
        # reference genome files
        File ref_dict
        File ref_fasta 
        File ref_fasta_index
        File bait_file
        File extra_amplicons_file
        File duplicate_amplicons_file
        # variant calling inputs
        File target_intervals
        File dbsnp
        File dbsnp_vcf_index
    }

    call FileUtils.FetchFilesTask as FetchCram {
            input:
                data_location = sample_data_location,
                file_types = ["cram"],
                recursive = false,
                file_match_keys = [subject_id, sample_id],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
    }

    call HaplotypeCallerGvcfGATK4.CramToBamTask {
            input:
                input_cram = select_first([FetchCram.cram]),
                sample_name = basename(select_first([FetchCram.cram]), ".cram"),
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
                samtools_path = "samtools"
    }

    call HaplotypeCallerGvcfGATK4.VariantCallingGenomePanelsTask {
        input:
            input_bam = CramToBamTask.output_bam,
            input_bam_index = CramToBamTask.output_bai,
            interval_list = target_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index
    }

    call LMMVariantReportTask {
        input:
            sample_id = sample_id,
            sample_run_id = sample_run_id,
            sample_lmm_id = sample_lmm_id,
            sample_barcode = sample_barcode,
            batch_id = batch_id,
            target_intervals = target_intervals,
            bait_file = bait_file,
            extra_amplicons_file = extra_amplicons_file,
            duplicate_amplicons_file = duplicate_amplicons_file,
            test_code = test_code,
            all_calls_vcf_file = VariantCallingGenomePanelsTask.all_calls_vcf_file,
            all_bases_noChr_vcf_file = VariantCallingGenomePanelsTask.all_bases_noChr_vcf_file,
            mgbpmbiofx_docker_image = orchutils_docker_image
    }

    output {

    }
}


task LMMVariantReportTask {
    input {
        String sample_id
        String sample_run_id
        String sample_lmm_id
        String sample_barcode
        String batch_id
        String test_code
        File target_intervals
        File bait_file
        File extra_amplicons_file
        File duplicate_amplicons_file
        File all_bases_noChr_vcf_file
        File all_calls_vcf_file
        String mgbpmbiofx_docker_image
        Int disk_size = 10
    }

    String lmm_vc_prefix = sample_run_id + '__' + sample_lmm_id + '__' + sample_id + '__' + sample_barcode + '__' + batch_id
    String xls_report_out = lmm_vc_prefix + '.variantReport.xls'
    String snps_out = lmm_vc_prefix + '.cleaned.common_ysnps.txt'
    String xml_report_out = lmm_vc_prefix + '.variantReport.xml'

    command <<<
        set -euxo pipefail
        $MGBPMBIOFXPATH/genome-panels/bin/lmm_variant_detection_report.py \
        ~{target_intervals} \
        ~{all_bases_noChr_vcf_file} \
        ~{all_calls_vcf_file} \
        --genome-build 'GRCh38' \
        --test-code ~{test_code} \
        --bait-set ~{bait_file} \
        --filter-on-test-code \
        --filter-on-gi-reference \
        '--webservice-url ' + config['ws_url'] \
        '--webservice-login ' + config['ws_login'] \
        '--webservice-pass ' + config['ws_pass'] \
        --run-id ~{sample_run_id} \
        --container-id ~{sample_lmm_id} \
        --batch-id ~{batch_id} \
        --accession-id ~{sample_id} \
        --index-barcode ~{sample_barcode} \
        --extra_amplicons_file ~{extra_amplicons_file} \
        --duplicate_amplicons_file ~{duplicate_amplicons_file} \
        --out ~{xls_report_out} \
        --snp ~{snps_out} \
        --xml ~{xml_report_out}
        
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        
    }

    output {

    }


}