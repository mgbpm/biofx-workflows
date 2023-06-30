version 1.0

import "./GATKWorkflow.wdl" as gatk

workflow RiskAllelesWorkflow {

    input {
        File input_cram
        File input_crai
        String sample_id
        String accession_id
        String test_code
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        File workflow_fileset
        String gcp_project_id
        String workspace_name
        String mgbpmbiofx_docker_image
    }

    String out_path = "outputs"
    
    call gatk.GATKWorkflow {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            sample_id =sample_id,
            test_code = test_code,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }


    call RiskTask {
        input:
            sample_id = sample_id,
            accession_id = accession_id,
            test_code = test_code,
            roi_bed = roi_bed,
            workflow_fileset = workflow_fileset,
            all_bases_vcf_file = GATKWorkflow.all_bases_vcf_file,
            all_bases_vcf_idx_file = GATKWorkflow.all_bases_vcf_idx_file,
            out_path = out_path,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    output {
        File genotype_xlsx = RiskTask.genotype_xlsx
        File genotype_txt = RiskTask.genotype_txt
        File risk_report = RiskTask.risk_report
    }
}

task RiskTask {
    input {
        String sample_id
        String accession_id
        String test_code
        File roi_bed
        File workflow_fileset
        File all_bases_vcf_file
        File all_bases_vcf_idx_file
        String out_path
        String gcp_project_id
        String workspace_name
        String mgbpmbiofx_docker_image
        Int disk_size = ceil(size(all_bases_vcf_file, "GB") + size(roi_bed, "GB")) + 10
    }

    command <<<
        set -euxo pipefail

        LIB_DIR="~{test_code}"_lib
        mkdir -p ~{out_path}/$LIB_DIR

        tar -xf "~{workflow_fileset}" -C ~{out_path}/$LIB_DIR
        
        ln -s ~{all_bases_vcf_file} ~{out_path}
        ln -s ~{all_bases_vcf_idx_file} ~{out_path}

        # Setup OMS client config
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/get-client-config.sh \
            -p ~{gcp_project_id} -w ~{workspace_name} -n oms > oms-client-config.json

        $MGBPMBIOFXPATH/biofx-risk-alleles/bin/run_risk.py \
            -s "~{sample_id}" \
            -a "~{accession_id}" \
            -t "~{test_code}" \
            -o "~{out_path}" \
            -roi "~{roi_bed}" \
            -l "~{out_path}/$LIB_DIR" \
            -oms oms-client-config.json
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk ~{disk_size} SSD"
    }

    output {
        File risk_report = "~{out_path}/~{sample_id}_~{test_code}.report.xlsx"
        File genotype_xlsx = "~{out_path}/~{sample_id}_~{test_code}.genotype.xlsx"
        File genotype_txt = "~{out_path}/~{sample_id}_~{test_code}.genotype.txt"
    }
}