version 1.0

import "./GATKWorkflow.wdl" as gatk

workflow PGxWorkflow {	

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

    call PGxTask {
        input:
            sample_id = sample_id,
            accession_id =accession_id,
            test_code = test_code,
            roi_bed = roi_bed,
            workflow_fileset = workflow_fileset,
            all_bases_vcf_file = GATKWorkflow.all_bases_vcf_file,
            all_bases_vcf_idx_file = GATKWorkflow.all_bases_vcf_idx_file,
            out_path = out_path,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    output {
        File FDA_report = PGxTask.FDA_report
        File CPIC_report = PGxTask.CPIC_report
        File genotype_xlsx = PGxTask.genotype_xlsx
        File genotype_txt = PGxTask.genotype_txt
    }
}


task PGxTask {
    input {
        String sample_id
        String accession_id
        String test_code
        File roi_bed
        File workflow_fileset
        File all_bases_vcf_file
        File all_bases_vcf_idx_file
        String out_path
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
        
        $MGBPMBIOFXPATH/biofx-pgx/bin/run_pgx.py \
            -s "~{sample_id}" \
            -a "~{accession_id}" \
            -t "~{test_code}" \
            -o "~{out_path}" \
            -roi "~{roi_bed}" \
            -l "~{out_path}/$LIB_DIR"
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        bootDiskSizeGb: 15
    }

    output {
        File CPIC_report = "~{out_path}/~{sample_id}_~{test_code}.CPIC_report.xlsx"
        File FDA_report = "~{out_path}/~{sample_id}_~{test_code}.FDA_report.xlsx"
        File genotype_xlsx = "~{out_path}/~{sample_id}_~{test_code}.genotype.xlsx"
        File genotype_txt = "~{out_path}/~{sample_id}_~{test_code}.genotype.txt"
    }
}
