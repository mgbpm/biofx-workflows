version 1.0

import "../../steps/FileUtils.wdl"

workflow SplitVcfFileToLocationWorkflow {
    input {
        String joint_vcf_file
        String target_location
        String mgbpmbiofx_docker_image
        String gcp_project_id
        String workspace_name
    }

    call FileUtils.FetchFilesTask {
        input:
            data_location = joint_vcf_file,
            file_types = ['vcf'],
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name
    }

    call SplitVcfFilesTask {
        input:
            joint_vcf_files = FetchFilesTask.all_files,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
    }

    scatter (vcf_file in SplitVcfFilesTask.individual_vcf_files) {
        String file_basename = basename(vcf_file)
        String sample_id = sub(file_basename, '\\.g?vcf[.bgz]*$', '')
        String vcf_file_str = vcf_file
        call FileUtils.CopyFilesTask {
            input:
                source_location = vcf_file_str,
                target_location = target_location + '/' + sample_id,
                flatten = true,
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name
        }
    }

    output {
        Array[String] target_files = flatten(CopyFilesTask.target_files)
    }
}

task SplitVcfFilesTask {
    input {
        Array[File] joint_vcf_files
        String mgbpmbiofx_docker_image
    }

    command <<<
        set -euxo pipefail

        mkdir bcf_files vcf_files
        for file in ~{sep=' ' joint_vcf_files}
        do
            bcftools +split $file -Ou -o bcf_files
        done

        for file in $(ls bcf_files)
        do
            bcftools view bcf_files/$file -Oz > vcf_files/$file
        done
        ls -1 vcf_files | sed 's#\(.*\)#vcf_files/\1#' > vcf_files.txt
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "4GB"
    }

    output {
        Array[File] individual_vcf_files = read_lines("vcf_files.txt")
    }
}