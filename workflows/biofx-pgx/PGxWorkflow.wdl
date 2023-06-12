version 1.0

import "./GATKWorkflow.wdl"

workflow PGxWorkflow {	

    input {
        File input_cram
        File input_crai
        String sample_id
        String accession_id
        String test_code
        String java_path = "/usr/lib/jvm/java-8-openjdk-amd64/bin/java"
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        File workflow_fileset
        String gatk_path = "/gatk/gatk"
        String mgbpmbiofx_docker_image
    }

    String out_path = "outputs"
    String out_prefix = sample_id + "_" + test_code
    String gvcf = out_path + '/' + out_prefix + ".g.vcf"
    String all_calls_vcf = out_path + '/' + out_prefix + '.allcalls.vcf'
    String ref_positions_vcf = out_path + '/' + out_prefix + '.ref_positions.vcf'
    String all_bases_vcf = out_path + '/' + out_prefix + '.allbases.vcf'

    call GATKWorkflow.HaplotypeCallerTask {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gatk_path = gatk_path,
            out_path = out_path,
            gvcf = gvcf,
            all_calls_vcf = all_calls_vcf
    }

    call GATKWorkflow.CreateRefSitesVCFTask {
        input:
            gvcf_file = HaplotypeCallerTask.gvcf_file,
            gvcf_idx_file = HaplotypeCallerTask.gvcf_idx_file,
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            all_calls_vcf_idx_file = HaplotypeCallerTask.all_calls_vcf_idx_file,
            ref_positions_vcf = ref_positions_vcf,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            out_path = out_path
    }

    call GATKWorkflow.SortVCFTask {
        input:
            gatk_path = gatk_path,
            java_path = java_path,
            out_path = out_path,
            ref_positions_vcf_file = CreateRefSitesVCFTask.ref_positions_vcf_file,
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            all_bases_vcf = all_bases_vcf,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
    }

    call PGxTask {
        input:
            sample_id = sample_id,
            accession_id =accession_id,
            test_code = test_code,
            roi_bed = roi_bed,
            workflow_fileset = workflow_fileset,
            gvcf_file = HaplotypeCallerTask.gvcf_file,
            gvcf_idx_file = HaplotypeCallerTask.gvcf_idx_file,
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            all_calls_vcf_idx_file = HaplotypeCallerTask.all_calls_vcf_idx_file,
            ref_positions_vcf_file = CreateRefSitesVCFTask.ref_positions_vcf_file,
            all_bases_vcf_file = SortVCFTask.all_bases_vcf_file,
            all_bases_vcf_idx_file = SortVCFTask.all_bases_vcf_idx_file,
            out_path = out_path,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    output {
        File all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file
        File all_calls_vcf_idx_file = HaplotypeCallerTask.all_calls_vcf_idx_file
        File gvcf_file = HaplotypeCallerTask.gvcf_file
        File gvcf_idx_file = HaplotypeCallerTask.gvcf_idx_file
        File ref_positions_vcf_file = CreateRefSitesVCFTask.ref_positions_vcf_file
        File all_bases_vcf_file = SortVCFTask.all_bases_vcf_file
        File all_bases_vcf_idx_file = SortVCFTask.all_bases_vcf_idx_file
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

        python3 $MGBPMBIOFXPATH/biofx-pgx/src/pgx.py \
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
        File CPIC_report = "~{out_path}/~{sample_id}_~{test_code}.CPIC_report.xlsx"
        File FDA_report = "~{out_path}/~{sample_id}_~{test_code}.FDA_report.xlsx"
        File genotype_xlsx = "~{out_path}/~{sample_id}_~{test_code}.genotype.xlsx"
        File genotype_txt = "~{out_path}/~{sample_id}_~{test_code}.genotype.txt"
    }
}
