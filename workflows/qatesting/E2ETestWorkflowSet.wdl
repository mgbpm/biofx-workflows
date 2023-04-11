version 1.0

import "../../steps/FileUtils.wdl"

workflow E2ETestWorkflowSet {
    input {
        Array[String] subject_id
        Array[String] sample_id
        Array[String] sample_data_location
        String mgbpmbiofx_docker_image
        String gcp_project_id
        String workspace_name
    }

    Array[Pair[String, String]] subject_sample_pairs = zip(subject_id, sample_id)
    Array[Pair[Pair[String, String], String]] subject_sample_dataloc_pairs = zip(subject_sample_pairs, sample_data_location)

    scatter (subject_sample_dataloc_pair in subject_sample_dataloc_pairs) {
        call FileUtils.FetchFilesTask {
            input:
                data_location = subject_sample_dataloc_pair.right,
                file_types = ['bam', 'cram'],
                docker_image = mgbpmbiofx_docker_image,
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
                subject_id = subject_sample_dataloc_pair.left.left,
                sample_id = subject_sample_dataloc_pair.left.right,
                bam = final_bam,
                bai = final_bai,
                cram = final_cram,
                crai = final_crai
        }
    }

    output {
        Array[File] vcf = MockVariantCallingAndCoverageTask.vcf
        Array[File] vcf_index = MockVariantCallingAndCoverageTask.vcf_index
        Array[File] cov_summary = MockVariantCallingAndCoverageTask.cov_summary
        Array[File] cov_gene_summary = MockVariantCallingAndCoverageTask.cov_gene_summary
        Array[File] cov_gene_details = flatten(MockVariantCallingAndCoverageTask.cov_gene_details)
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