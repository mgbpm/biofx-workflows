version 1.0

workflow TestWorkflowFailure {
    input {
        String subject_id
    }

    call FailureTask {
        input:
            subject_id = subject_id
    }

    output {
        String output_subject_id = FailureTask.output_subject_id
    }
}

task FailureTask {
    input {
        String subject_id
    }

    command <<<
        set -euxo pipefail

        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: "4GB"
    }

    output {
        String output_subject_id = "~{subject_id}"
    }
}