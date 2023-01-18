version 1.0

workflow UnitTestWorkflowMany {
    input {
        Array[String?] input_1
        Array[String?] input_2
        Array[String?] input_3
        Array[String?] input_4
        Array[String?] input_5
    }

    call UnitTestWorkflowManyTask {
        input:
            input_1 = input_1,
            input_2 = input_2,
            input_3 = input_3,
            input_4 = input_4,
            input_5 = input_5,
    }

    output {
        String output_1 = UnitTestWorkflowManyTask.output_1
        String output_2 = UnitTestWorkflowManyTask.output_2
        String output_3 = UnitTestWorkflowManyTask.output_3
        String output_4 = UnitTestWorkflowManyTask.output_4
        String output_5 = UnitTestWorkflowManyTask.output_5
        File outputs_file = UnitTestWorkflowManyTask.outputs_file
    }
}

task UnitTestWorkflowManyTask {
    input {
        Array[String?] input_1
        Array[String?] input_2
        Array[String?] input_3
        Array[String?] input_4
        Array[String?] input_5
    }

    command <<<
        set -euxo pipefail

        echo "input_1=~{sep=' ' input_1}" >> outputs.txt
        echo "input_2=~{sep=' ' input_2}" >> outputs.txt
        echo "input_3=~{sep=' ' input_3}" >> outputs.txt
        echo "input_4=~{sep=' ' input_4}" >> outputs.txt
        echo "input_5=~{sep=' ' input_5}" >> outputs.txt

        INPUT1="~{sep=' ' input_1}"
        INPUT2="~{sep=' ' input_2}"
        INPUT3="~{sep=' ' input_3}"
        INPUT4="~{sep=' ' input_4}"
        INPUT5="~{sep=' ' input_5}"

        [ -z "${INPUT1}" ] && INPUT1="Empty"
        [ -z "${INPUT2}" ] && INPUT2="Empty"
        [ -z "${INPUT3}" ] && INPUT3="Empty"
        [ -z "${INPUT4}" ] && INPUT4="Empty"
        [ -z "${INPUT5}" ] && INPUT5="Empty"

        echo "${INPUT1}" >> input_1.txt
        echo "${INPUT2}" >> input_2.txt
        echo "${INPUT3}" >> input_3.txt
        echo "${INPUT4}" >> input_4.txt
        echo "${INPUT5}" >> input_5.txt
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: "4GB"
    }

    output {
        String output_1 = read_string("input_1.txt")
        String output_2 = read_string("input_2.txt")
        String output_3 = read_string("input_3.txt")
        String output_4 = read_string("input_4.txt")
        String output_5 = read_string("input_5.txt")
        File outputs_file = "outputs.txt"
    }
}