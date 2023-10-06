version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/HaplotypeCallerGvcfGATK4.wdl"
import "../../steps/IgvReport.wdl"
import "../../steps/Utilities.wdl"

workflow BgwgsPanelWorkflow {

    input {
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id
        String workspace_name
        # Docker images
        String orchutils_docker_image
        String genome_panels_docker_image
        # Run, sample id and data location
        String subject_id
        String sample_id
        String sample_lmm_id
        String sample_barcode
        String batch_id
        String sample_data_location
        Boolean fetch_cram = true
        Boolean fetch_bam = true
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
        # Reporting steps
        String igvreport_docker_image
    }

    # Prefer a BAM file to avoid conversion
    if (fetch_bam) {
        call FileUtils.FetchFilesTask as FetchBam {
            input:
                data_location = sample_data_location,
                file_types = ["bam"],
                recursive = false,
                file_match_keys = [sample_lmm_id, subject_id],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
        }
    }

    # no BAM, fallback to CRAM
    if (fetch_cram && (!defined(FetchBam.bam) || !defined(FetchBam.bai))) {
        call FileUtils.FetchFilesTask as FetchCram {
            input:
                data_location = sample_data_location,
                file_types = ["cram"],
                recursive = false,
                file_match_keys = [sample_lmm_id, subject_id],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                disk_size = fetch_disk_size
        }
    }

    # Validate that we got workable files from the data location
    if (!defined(FetchBam.bam) && !defined(FetchCram.cram)) {
        call Utilities.FailTask as MissingBamOrCramFailure {
            input:
                error_message = "BAM or CRAM file not found in ~{sample_data_location}"
        }
    }
    if (defined(FetchCram.cram) && !defined(FetchCram.crai)) {
        call Utilities.FailTask as MissingCraiFailure {
            input:
                error_message = "Index file for CRAM " + basename(select_first([FetchCram.cram])) + " not found"
        }
    }
    if (!defined(FetchCram.cram) && defined(FetchBam.bam) && !defined(FetchBam.bai)) {
        call Utilities.FailTask as MissingBaiFailure {
            input:
                error_message = "Index file for BAM " + basename(select_first([FetchCram.bam])) + " not found"
        }
    }
    
    # Convert CRAM to BAM
    if (defined(FetchCram.cram) && defined(FetchCram.crai)) {
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
    }

    # final inputs - either CRAM or BAM
    File sample_bam = select_first([FetchBam.bam, CramToBamTask.output_bam])
    File sample_bai = select_first([FetchBam.bai, CramToBamTask.output_bai])

    # Run haplotype caller
    call HaplotypeCallerGvcfGATK4.GenomePanelsVariantCallingTask {
        input:
            input_bam = sample_bam,
            input_bam_index = sample_bai,
            interval_list = target_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index
    }

    call CreateRefSitesTask {
        input:
            gvcf_file = GenomePanelsVariantCallingTask.gvcf_file,
            all_calls_vcf_file = GenomePanelsVariantCallingTask.all_calls_vcf_file,
            docker = genome_panels_docker_image 
    }

    call HaplotypeCallerGvcfGATK4.GenomePanelsRefSitesSortTask {
        input:
            gvcf_file = GenomePanelsVariantCallingTask.gvcf_file,
            all_calls_vcf_file = GenomePanelsVariantCallingTask.all_calls_vcf_file,
            ref_positions_vcf_file = CreateRefSitesTask.ref_positions_vcf_file
    }

    call LMMVariantReportTask {
        input:
            subject_id = subject_id,
            sample_id = sample_id,
            sample_lmm_id = sample_lmm_id,
            sample_barcode = sample_barcode,
            batch_id = batch_id,
            target_intervals = target_intervals,
            bait_file = bait_file,
            extra_amplicons_file = extra_amplicons_file,
            duplicate_amplicons_file = duplicate_amplicons_file,
            test_code = test_code,
            all_calls_vcf_file = GenomePanelsVariantCallingTask.all_calls_vcf_file,
            all_bases_vcf_file = GenomePanelsRefSitesSortTask.all_bases_vcf_file,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            docker_image = genome_panels_docker_image
    }

    call CreateBedRegionsFromXLSTask {
        input:
            lmm_xls_report_file = LMMVariantReportTask.xls_report_out_file,
            docker_image = genome_panels_docker_image
    }


    call IgvReport.IgvReportFromGenomePanelsBedTask {
        input:
            bed_file = CreateBedRegionsFromXLSTask.output_bed_file,
            sample_bam = sample_bam,
            sample_bai = sample_bai,
            output_basename = sample_id + '__' + sample_lmm_id + '__' + subject_id + '__' + sample_barcode + '__' + batch_id + ".igvreport",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            docker_image = igvreport_docker_image
    }

    output {
        # Variant calling outputs
        File all_calls_vcf_file = GenomePanelsVariantCallingTask.all_calls_vcf_file 
        File all_calls_vcf_idx_file = GenomePanelsVariantCallingTask.all_calls_vcf_idx_file
        File gvcf_file = GenomePanelsVariantCallingTask.gvcf_file 
        File gvcf_idx_file = GenomePanelsVariantCallingTask.gvcf_idx_file
        File ref_positions_vcf_file = CreateRefSitesTask.ref_positions_vcf_file
        File all_bases_vcf_file = GenomePanelsRefSitesSortTask.all_bases_vcf_file
        File all_bases_vcf_idx_file = GenomePanelsRefSitesSortTask.all_bases_vcf_idx_file

        # LMM variant detection outputs
        File all_bases_noChr_vcf_file = LMMVariantReportTask.all_bases_noChr_vcf_file
        File? snps_out_file = LMMVariantReportTask.snps_out_file
        File xls_report_out_file = LMMVariantReportTask.xls_report_out_file
        File xml_report_out_file = LMMVariantReportTask.xml_report_out_file

        # LMM IGV reports
        File follow_up_regions_file = CreateBedRegionsFromXLSTask.output_bed_file
        File lmm_igv_report = IgvReportFromGenomePanelsBedTask.igv_report
    }
}


task CreateRefSitesTask {
  input {
        # Command parameters
        File gvcf_file
        File all_calls_vcf_file
        String docker
        Int? mem_gb
        Int? disk_space_gb
        Boolean use_ssd = false
        Int? preemptible_attempts
    }

    Int machine_mem_gb = select_first([mem_gb, 24])

    Int disk_size = ceil(size(gvcf_file, "GB") + size(all_calls_vcf_file, "GB")) + 10

    String sample_basename = basename(gvcf_file, ".g.vcf")
    String ref_positions_vcf = sample_basename + ".ref_positions.vcf"

    command <<<
      set -euxo pipefail

      "$MGBPMBIOFXPATH/biofx-genome-panels/bin/create_ref_sites_vcf.py" \
      -g "~{gvcf_file}" \
      -c "~{all_calls_vcf_file}" \
      -o "~{ref_positions_vcf}"
    >>>

    runtime {
      docker: docker
      memory: machine_mem_gb + " GB"
      disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
      preemptible: select_first([preemptible_attempts, 2])
    }

    output {
      File ref_positions_vcf_file = "~{ref_positions_vcf}"
    }
}


task LMMVariantReportTask {
    input {
        String subject_id
        String sample_id
        String sample_lmm_id
        String sample_barcode
        String batch_id
        String test_code
        File target_intervals
        File bait_file
        File extra_amplicons_file
        File duplicate_amplicons_file
        File all_bases_vcf_file
        File all_calls_vcf_file
        String gcp_project_id
        String workspace_name
        String docker_image
        Int disk_size = 10
    }

    String sample_basename = basename(all_bases_vcf_file, ".allbases.vcf")
    String all_bases_noChr_vcf = sample_basename + ".allbases.noChr.vcf"
    String lmm_vc_prefix = sample_id + '__' + sample_lmm_id + '__' + subject_id + '__' + sample_barcode + '__' + batch_id
    String xls_report_out = lmm_vc_prefix + '.variantReport.xls'
    String snps_out = lmm_vc_prefix + '.cleaned.common_ysnps.txt'
    String xml_report_out = lmm_vc_prefix + '.variantReport.xml'

    command <<<
        set -euxo pipefail

        #make a copy of three files for lmm_variant_detection_report.py
        sed 's/chr//' ~{all_bases_vcf_file} >  ~{all_bases_noChr_vcf}

        "$MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh" \
            -p ~{gcp_project_id} -w ~{workspace_name} -n gil > gil-client-config.json

        "$MGBPMBIOFXPATH/biofx-genome-panels/bin/lmm_variant_detection_report.py" \
        ~{target_intervals} \
        ~{all_bases_noChr_vcf} \
        ~{all_calls_vcf_file} \
        --genome-build 'GRCh38' \
        --test-code ~{test_code} \
        --bait-set ~{bait_file} \
        --filter-on-test-code \
        --filter-on-gi-reference \
        --run-id ~{sample_id} \
        --container-id ~{sample_lmm_id} \
        --batch-id ~{batch_id} \
        --accession-id ~{subject_id} \
        --index-barcode ~{sample_barcode} \
        --extra_amplicons_file ~{extra_amplicons_file} \
        --duplicate_amplicons_file ~{duplicate_amplicons_file} \
        --out ~{xls_report_out} \
        --snp ~{snps_out} \
        --xml ~{xml_report_out} \
        --webservice-config gil-client-config.json 
        
    >>>

    runtime {
        docker: "~{docker_image}"
        # memory: machine_mem_gb + " GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 2
    }

    output {
        File all_bases_noChr_vcf_file = "~{all_bases_noChr_vcf}"
        File? snps_out_file = "~{snps_out}"
        File xls_report_out_file = "~{xls_report_out}"
        File xml_report_out_file = "~{xml_report_out}"
    }
}

task CreateBedRegionsFromXLSTask {
    input {
        File lmm_xls_report_file
        String output_bed = "follow_up_regions.bed"
        String docker_image
        Int disk_size = 10
    }

    command <<<
        set -euxo pipefail

        "$MGBPMBIOFXPATH/biofx-genome-panels/bin/igv_lmm_wrapper.py" \
        ~{lmm_xls_report_file} \
        ~{output_bed}
    >>>

    runtime {
        docker: "~{docker_image}"
        # memory: machine_mem_gb + " GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 2
    }

    output {
        File output_bed_file = "~{output_bed}"
    }
}