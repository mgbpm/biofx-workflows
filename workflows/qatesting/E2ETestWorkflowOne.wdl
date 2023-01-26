version 1.0

import "../../steps/FileUtils.wdl"

workflow E2ETestWorkflowOne {
    input {
        String subject_id
        String sample_id
        String sample_data_location
        String mgbpmbiofx_docker_image
        String gcp_project_id
        String workspace_name
    }

    call FileUtils.FetchFilesTask {
        input:
            data_location = sample_data_location,
            file_types = ['bam', 'cram'],
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name
    }

    # type coercion will cause failure if task outputs are not defined
    File final_bam = select_first([FetchFilesTask.bam])
    File final_bai = select_first([FetchFilesTask.bai])
    File final_cram = select_first([FetchFilesTask.cram])
    File final_crai = select_first([FetchFilesTask.crai])

    call MockVariantCallingAndCoverageTask {
        input:
            subject_id = subject_id,
            sample_id = sample_id,
            bam = final_bam,
            bai = final_bai,
            cram = final_cram,
            crai = final_crai
    }

    output {
        File vcf = MockVariantCallingAndCoverageTask.vcf
        File vcf_index = MockVariantCallingAndCoverageTask.vcf_index
        File cov_summary = MockVariantCallingAndCoverageTask.cov_summary
        File cov_gene_summary = MockVariantCallingAndCoverageTask.cov_gene_summary
        Array[File] cov_gene_details = MockVariantCallingAndCoverageTask.cov_gene_details
    }
}

task MockVariantCallingAndCoverageTask {
    input {
        String subject_id
        String sample_id
        File bam
        File bai
        File cram
        File crai
    }

    command <<<
        set -euxo pipefail

        # validate inputs
        [ -z "~{subject_id}" ] && exit 1
        [ -z "~{sample_id}" ] && exit 1
        [ ! -r "~{bam}" -a ! -r "~{cram}" ] && exit 1
        [ -r "~{bam}" -a ! -r "~{bai}" ] && exit 1
        [ -r "~{cram}" -a ! -r "~{crai}" ] && exit 1

        # create output files
        echo "THIS IS THE VCF FILE" > "~{subject_id}_~{sample_id}.vcf.gz"
        echo "THIS IS THE VCF INDEX FILE" > "~{subject_id}_~{sample_id}.vcf.gz.idx"
        echo "THIS IS THE COVERAGE SUMMARY FILE" > "~{subject_id}_~{sample_id}.cov.summary.txt"
        echo "THIS IS THE GENE COVERAGE SUMMARY FILE" > "~{subject_id}_~{sample_id}.cov.gene_summary.txt"
        echo "THIS IS THE COVERAGE DETAILS FOR ABC1" > "~{subject_id}_~{sample_id}.cov.ABC1.txt"
        echo "THIS IS THE COVERAGE DETAILS FOR ABC2" > "~{subject_id}_~{sample_id}.cov.ABC2.txt"
        echo "THIS IS THE COVERAGE DETAILS FOR ABC3" > "~{subject_id}_~{sample_id}.cov.ABC3.txt"
        echo "THIS IS THE COVERAGE DETAILS FOR ABC4" > "~{subject_id}_~{sample_id}.cov.ABC4.txt"
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: "4GB"
    }

    output {
        File vcf = "~{subject_id}_~{sample_id}.vcf.gz"
        File vcf_index = "~{subject_id}_~{sample_id}.vcf.gz.idx"
        File cov_summary = "~{subject_id}_~{sample_id}.cov.summary.txt"
        File cov_gene_summary = "~{subject_id}_~{sample_id}.cov.gene_summary.txt"
        Array[File] cov_gene_details = ["~{subject_id}_~{sample_id}.cov.ABC1.txt", "~{subject_id}_~{sample_id}.cov.ABC2.txt", "~{subject_id}_~{sample_id}.cov.ABC3.txt", "~{subject_id}_~{sample_id}.cov.ABC4.txt"]
    }
}