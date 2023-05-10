version 1.0

workflow PGxWorkflow {	

    input {
        File input_cram
        File input_crai
        String sample_id
        String accession_id
        String test_code
        String java_options = "-Xms12g -Xmx40g"
        String java_path = "/usr/lib/jvm/java-8-openjdk-amd64/bin/java"
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File workflow_fileset
        String gatk_path = "/gatk/gatk"
        String mgbpmbiofx_docker_image
    }

    call GATKTask {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            sample_id = sample_id,
            test_code = test_code,
            java_options = java_options,
            java_path = java_path,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gatk_path = gatk_path
    }

    call PGxTask {
        input:
            sample_id = sample_id,
            accession_id =accession_id,
            test_code = test_code,
            roi_bed = roi_bed,
            workflow_fileset = workflow_fileset,
            supporting_files = GATKTask.supporting_files,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    output {
        File FDA_report = PGxTask.FDA_report
        File CPIC_report = PGxTask.CPIC_report
        File genotype_xlsx = PGxTask.genotype_xlsx
    }
}



# Task for gatk haplotype caller
task GATKTask {
    input {
        File input_cram
        File input_crai
        String sample_id
        String test_code
        String java_options
        String java_path
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        String mgbpmbiofx_docker_image
        String gatk_path
    }

    command <<<
        set -euxo pipefail

        mkdir outputs

        python3 $MGBPMBIOFXPATH/biofx-pgx/src/genotype_variants.py \
        -b "~{input_cram}" \
        -s "~{sample_id}" \
        -t "~{test_code}" \
        -o outputs \
        -R "~{reference_fasta}" \
        -roi "~{roi_bed}" \
        -db "~{dbsnp}" \
        -jp "~{java_path}" \
        -jo "~{java_options}" \
        -gp "~{gatk_path}"
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk 100 SSD"
    }

    output {
        Array[File]+ supporting_files = glob("outputs/~{sample_id}/~{test_code}/supporting/*")
    }
}

task PGxTask {
    input {
        String sample_id
        String accession_id
        String test_code
        File roi_bed
        File workflow_fileset
        Array[File]+ supporting_files
        String mgbpmbiofx_docker_image
    }

    command <<<
        set -euxo pipefail
        LIB_DIR="~{test_code}"_lib
        mkdir -p outputs/"~{sample_id}"/"~{test_code}"/supporting/$LIB_DIR

        tar -xf "~{workflow_fileset}" -C outputs/"~{sample_id}"/"~{test_code}"/supporting/$LIB_DIR
        
        for x in ~{sep=' ' supporting_files}
        do
            cp "${x}" outputs/"~{sample_id}"/"~{test_code}"/supporting/
        done;
        
        python3 $MGBPMBIOFXPATH/biofx-pgx/src/pgx.py \
            -s "~{sample_id}" \
            -a "~{accession_id}" \
            -t "~{test_code}" \
            -o outputs \
            -roi "~{roi_bed}" \
            -l outputs/"~{sample_id}"/"~{test_code}"/supporting/$LIB_DIR
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk 100 SSD"
    }

    output {
        File CPIC_report = "outputs/~{sample_id}/~{test_code}/~{sample_id}_~{test_code}.CPIC_report.xlsx"
        File FDA_report = "outputs/~{sample_id}/~{test_code}/~{sample_id}_~{test_code}.FDA_report.xlsx"
        File genotype_xlsx = "outputs/~{sample_id}/~{test_code}/~{sample_id}_~{test_code}.genotype.xlsx"
    }
}
