version 1.0

workflow PGxWorkflow {	

    input {
        File input_cram
        String sample_id
        String accession_id
        String test_code
        String java_options = "-Xms12g -Xmx40g"
        String java_path = "/usr/lib/jvm/java-8-openjdk-amd64/bin/java"
        File reference_fasta
        File roi_bed
        File dbsnp
        File workflow_fileset
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
        String gatk_path = "/gatk/gatk"
        String mgbpmbiofx_docker_image
    }

    call GATKTask {
        input:
            input_cram = input_cram,
            sample_id = sample_id,
            test_code = test_code,
            java_options = java_options,
            java_path = java_path,
            reference_fasta = reference_fasta,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            gatk_docker = gatk_docker,
            gatk_path = gatk_path
    }

    call PGxTask {
        input:
            input_cram = input_cram,
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
        String sample_id
        String test_code
        String java_options
        String java_path
        File reference_fasta
        File roi_bed
        File dbsnp
        String gatk_docker
        String gatk_path
    }

    command <<<
        set -euxo pipefail
        
        python3 $MGBPMBIOFXPATH/biofx-pgx/src/genotype_variants.py \
        -b ~{input_cram} \
        -s ~{sample_id} \
        -t ~{test_code} \
        -R ~{reference_fasta} \
        -roi ~{roi_bed} \
        -db ~{dbsnp} \
        -jp ~{java_path} \
        -jo ~{java_options}\
        -gp ~{gatk_path}
    >>>

    runtime {
        docker: ~{gatk_docker}
        # Check log file to see memory requirements 
        # memory: "4GB"
    }

    output {
        Array[File]+ supporting_files = glob("~{sample_id}/~{test_code}/supporting/*")
    }
}

task PGxTask {
    input {
        File input_cram
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
        LIB_DIR=$MGBPMBIOFXPATH/"~{test_code}"_lib
        mkdir -p $LIB_DIR

        tar -xf "~{workflow_fileset}" -C $LIB_DIR
        
        $MGBPMBIOFXPATH/biofx-pgx/src/pgx.py \
            -b "~{input_cram}" \
            -s "~{sample_id}" \
            -a "~{accession_id}" \
            -t "~{test_code}" \
            -roi "~{roi_bed}" \
            -l $LIB_DIR
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "4GB"
    }

    output {
        File CPIC_report = "~{sample_id}/~{test_code}/~{sample_id}_~{test_code}.CPIC_report.xlsx"
        File FDA_report = "~{sample_id}/~{test_code}/~{sample_id}_~{test_code}.FDA_report.xlsx"
        File genotype_xlsx = "~{sample_id}/~{test_code}/~{sample_id}_~{test_code}.genotype.xlsx"
        Array[File]+ supporting_files = glob("~{sample_id}/~{test_code}/supporting/*")
    }
}
