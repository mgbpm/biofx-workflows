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

task MapSampleTrackerSexToPlinkSexTask {
    input {
        Array[String] sample_tracker_sexes
    }

    command <<<
        set -euxo pipefail

        for sex in ~{sep=' ' sample_tracker_sexes}
        do
            plink_sex="0"
            case "$sex" in
                "Unspecified" | "unspecified" | "Intersex" | "intersex" | "368000002" | "368000003" | "0" | "U")
                    plink_sex="0"
                    ;;
                "Male" | "male" | "MALE" | "368000001" | "1" | "M")
                    plink_sex="1"
                    ;;
                "Female" | "female" | "FEMALE" | "368000000" | "F" | "2")
                    plink_sex="2"
                    ;;
            esac
            echo "$plink_sex" >> output.txt
        done
    >>>


    runtime {
        docker: "ubuntu:latest"
    }

    output {
        Array[String] plink_sexes = read_lines("output.txt")
    }
}