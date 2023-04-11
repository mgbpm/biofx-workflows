version 1.0

task FailTask {
    input {
        String error_message = "workflow failed"
    }

    command <<<
        echo "ERROR: ~{error_message}" 1>&2
        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
    }

    output {

    }
}