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
            gvcf_file = GATKWorkflow.gvcf_file,
            gvcf_idx_file = GATKWorkflow.gvcf_idx_file,
            all_calls_vcf_file = GATKWorkflow.all_calls_vcf_file,
            all_calls_vcf_idx_file = GATKWorkflow.all_calls_vcf_idx_file,
            ref_positions_vcf_file = GATKWorkflow.ref_positions_vcf_file,
            all_bases_vcf_file = GATKWorkflow.all_bases_vcf_file,
            all_bases_vcf_idx_file = GATKWorkflow.all_bases_vcf_idx_file,
            out_path = out_path,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    output {
        File all_calls_vcf_file = GATKWorkflow.all_calls_vcf_file
        File all_calls_vcf_idx_file = GATKWorkflow.all_calls_vcf_idx_file
        File gvcf_file = GATKWorkflow.gvcf_file
        File gvcf_idx_file = GATKWorkflow.gvcf_idx_file
        File ref_positions_vcf_file = GATKWorkflow.ref_positions_vcf_file
        File all_bases_vcf_file = GATKWorkflow.all_bases_vcf_file
        File all_bases_vcf_idx_file = GATKWorkflow.all_bases_vcf_idx_file
        File FDA_report = RiskTask.risk_report
    }
}

task RiskTask {
    input {
        String sample_id
        String accession_id
        String test_code
        File roi_bed
        File workflow_fileset
        File gvcf_file
        File gvcf_idx_file
        File all_calls_vcf_file
        File all_calls_vcf_idx_file
        File ref_positions_vcf_file
        File all_bases_vcf_file
        File all_bases_vcf_idx_file
        String out_path
        String mgbpmbiofx_docker_image
    }

    command <<<
        set -euxo pipefail

        LIB_DIR="~{test_code}"_lib
        mkdir -p ~{out_path}/$LIB_DIR

        tar -xf "~{workflow_fileset}" -C ~{out_path}/$LIB_DIR
        
        ln -s ~{all_calls_vcf_file} ~{out_path}
        ln -s ~{all_calls_vcf_idx_file} ~{out_path}
        ln -s ~{all_bases_vcf_file} ~{out_path}
        ln -s ~{all_bases_vcf_idx_file} ~{out_path}
        ln -s ~{gvcf_file} ~{out_path}
        ln -s ~{gvcf_idx_file} ~{out_path}
        ln -s ~{ref_positions_vcf_file} ~{out_path}

        python $MGBPMBIOFXPATH/biofx-risk-alleles/src/risk.py \
            -s "~{sample_id}" \
            -a "~{accession_id}" \
            -t "~{test_code}" \
            -o "~{out_path}" \
            -roi "~{roi_bed}" \
            -l "~{out_path}/$LIB_DIR"
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk 100 SSD"
    }

    output {
        File risk_report = "~{out_path}/~{sample_id}_~{test_code}.report.xlsx"
        File genotype_xlsx = "~{out_path}/~{sample_id}_~{test_code}.genotype.xlsx"
        File genotype_txt = "~{out_path}/~{sample_id}_~{test_code}.genotype.txt"
    }
}