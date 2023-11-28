version 1.0

import "../../steps/QCEval.wdl" # for qc of each vcf
import "../../steps/FASTUtils.wdl" # for loading to FAST
import "../../steps/Utilities.wdl" # for fail task

workflow SampleLoadingWorkflow {
    input {
        # Dataset and sample inputs
        File input_vcf
        String sample_ID
        String dataset
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Orchestration utils docker image
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230921"
        # reference genome
        String reference_build = "GRCh38"
        # qceval inputs
        String qceval_project_type
        String qceval_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005"
        # FAST loading inputs
        String queue_folder = "gs://biofx-task-queue/biobank-ror-1"
        Int queue_limit = 200
        Int queue_wait_limit_hrs = 15
        Int fast_data_load_wait_interval_secs = 120
        Int fast_data_load_wait_max_intervals = 300
        Boolean has_haploid_sites = false
        String sample_data_load_config_name = "Sample_VCF_PPM_Eval"
    }

    # Annotate individual sample VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = input_vcf,
            project_type = qceval_project_type,
            output_basename = dataset + "_" + sample_ID + ".qceval",
            docker_image = qceval_docker_image
    }
    # Enter the queue for loading sample data
    call EnterQueueTask {
        input:
            input_vcf = QCEvalTask.output_vcf_gz,
            gcp_project_id = gcp_project_id,
            queue_folder = queue_folder,
            queue_limit = queue_limit,
            queue_wait_limit_hrs = queue_wait_limit_hrs,
            docker_image = orchutils_docker_image
    }
    # Load sample data to FAST
    if (defined(EnterQueueTask.queue_entry_id)) {
        call FASTUtils.FASTDataLoadTask as QCEvalLoadTask {
            input:
                reference_build = reference_build,
                vcf_file = QCEvalTask.output_vcf_gz,
                has_haploid_sites = has_haploid_sites,
                sample_data_name = dataset + "_" + sample_ID,
                data_load_config_name = sample_data_load_config_name,
                data_load_target = "SAMPLE_DATA",
                annotation_record_ts = "now",
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
        # Wait for data loads to complete
        call FASTUtils.FASTWaitForDataLoadsTask {
            input:
                data_load_ids = [ QCEvalLoadTask.data_load_id ],
                check_interval_secs = fast_data_load_wait_interval_secs,
                max_checks = fast_data_load_wait_max_intervals,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name,
                docker_image = orchutils_docker_image
        }
        if (defined(FASTWaitForDataLoadsTask.wait_result.total)) {
            # Leave queue when data load is successful
            call LeaveQueueTask {
                input:
                    gcp_project_id = gcp_project_id,
                    queue_folder = queue_folder,
                    queue_entry_id = EnterQueueTask.queue_entry_id,
                    docker_image = orchutils_docker_image
            }
        }
        if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.completed) {
            call Utilities.FailTask as DataLoadTimeout {
                input:
                    error_message = "One or more FAST data loads did not complete within allotted time"
            }
        }
        if (FASTWaitForDataLoadsTask.wait_result.total != FASTWaitForDataLoadsTask.wait_result.succeeded) {
            call Utilities.FailTask as DataLoadFailure {
                input:
                    error_message = "One or more FAST data loads did not succeed"
            }
        }
        
    }

    output {
        # Annotated qceval VCF
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        # FAST sample data name
        String fast_sample_data_name = dataset + "_" + sample_ID
        String queue_entry_id = EnterQueueTask.queue_entry_id
    }
}

task EnterQueueTask {
    input {
        File input_vcf
        String gcp_project_id
        String queue_folder
        Int queue_limit
        Int queue_wait_limit_hrs
        String docker_image
        Int disk_size = 5 + ceil(size(input_vcf, "GB"))
    }

    command <<<
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w global -r "~{queue_folder}"
        
        # Enter the task queue
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/enter_task_queue.py --queue-folder "~{queue_folder}" \
            --entry-details "source = BiobankRoRSampleLoadingWorkflow.wdl, vcf = ~{input_vcf}" \
            --queue-limit ~{queue_limit} \
            --wait --wait-ttl-hrs ~{queue_wait_limit_hrs} --wait-interval-mins 5 | tee queue_entry_id
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
    }

    output {
        String queue_entry_id = read_string("queue_entry_id")
    }
}

task LeaveQueueTask {
    input {
        String gcp_project_id
        String queue_folder
        String queue_entry_id
        String docker_image
        Int disk_size = 5
        Int preemptible = 1
    }

    command <<<
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p "~{gcp_project_id}" -w global -r "~{queue_folder}"

        # Remove the queue entry
        $MGBPMBIOFXPATH/biofx-orchestration-utils/bin/remove_task_queue_entry.py --queue-folder "~{queue_folder}" --entry-id "~{queue_entry_id}"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: preemptible
    }

    output {
        # Empty output
    }
}