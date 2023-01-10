version 1.0

workflow UnitTestWorkflowOne {
    input {
        String input_1
        String input_2
        String input_3
        String input_4
        String input_5
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
        String input_1
        String input_2
        String input_3
        String input_4
        String input_5
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
        String output_1 = "input_1=" + input_1
        String output_2 = "input_2=" + input_2
        String output_3 = "input_3=" + input_3
        String output_4 = "input_4=" + input_4
        String output_5 = "input_5=" + input_5
        File outputs_file = "outputs.txt"
    }
}