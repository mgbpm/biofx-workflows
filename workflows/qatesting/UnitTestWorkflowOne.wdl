version 1.0

workflow UnitTestWorkflowOne {
    input {
        String? input_1
        String? input_2
        String? input_3
        String? input_4
        String? input_5
    }

    call UnitTestWorkflowOneTask {
        input:
            input_1 = input_1,
            input_2 = input_2,
            input_3 = input_3,
            input_4 = input_4,
            input_5 = input_5,
    }

    output {
        String output_1 = UnitTestWorkflowOneTask.output_1
        String output_2 = UnitTestWorkflowOneTask.output_2
        String output_3 = UnitTestWorkflowOneTask.output_3
        String output_4 = UnitTestWorkflowOneTask.output_4
        String output_5 = UnitTestWorkflowOneTask.output_5
        File outputs_file = UnitTestWorkflowOneTask.outputs_file
    }
}

task UnitTestWorkflowOneTask {
    input {
        String? input_1
        String? input_2
        String? input_3
        String? input_4
        String? input_5
    }

    command <<<
        set -euxo pipefail

        echo "input_1=~{input_1}" >> outputs.txt
        echo "input_2=~{input_2}" >> outputs.txt
        echo "input_3=~{input_3}" >> outputs.txt
        echo "input_4=~{input_4}" >> outputs.txt
        echo "input_5=~{input_5}" >> outputs.txt
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: "4GB"
    }

    output {
        String output_1 = if defined(input_1) then "input_1=" + select_first([input_1]) else "input_1=Empty"
        String output_2 = if defined(input_1) then "input_2=" + select_first([input_2]) else "input_2=Empty"
        String output_3 = if defined(input_1) then "input_3=" + select_first([input_3]) else "input_3=Empty"
        String output_4 = if defined(input_1) then "input_4=" + select_first([input_4]) else "input_4=Empty"
        String output_5 = if defined(input_1) then "input_5=" + select_first([input_5]) else "input_5=Empty"
        File outputs_file = "outputs.txt"
    }
}