version 1.0

workflow UnitTestWorkflowTwo {
    input {
        Array[String?] input_1
        Array[String?] input_2
        Array[String?] input_3
        Array[String?] input_4
        Array[String?] input_5
    }

    call UnitTestWorkflowTwoTask {
        input:
            input_1 = input_1,
            input_2 = input_2,
            input_3 = input_3,
            input_4 = input_4,
            input_5 = input_5,
    }

    output {
        Array[String] output_1 = UnitTestWorkflowTwoTask.output_1
        Array[String] output_2 = UnitTestWorkflowTwoTask.output_2
        Array[String] output_3 = UnitTestWorkflowTwoTask.output_3
        Array[String] output_4 = UnitTestWorkflowTwoTask.output_4
        Array[String] output_5 = UnitTestWorkflowTwoTask.output_5
        File outputs_file = UnitTestWorkflowTwoTask.outputs_file
    }
}

task UnitTestWorkflowTwoTask {
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
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: "4GB"
    }

    output {
        Array[String] output_1 = if defined(input_1) then ["input_1=" + select_all(input_1)] else ["input_1=Empty"]
        Array[String] output_2 = if defined(input_2) then ["input_2=" + select_all(input_2)] else ["input_2=Empty"]
        Array[String] output_3 = if defined(input_3) then ["input_3=" + select_all(input_3)] else ["input_3=Empty"]
        Array[String] output_4 = if defined(input_4) then ["input_4=" + select_all(input_4)] else ["input_4=Empty"]
        Array[String] output_5 = if defined(input_5) then ["input_5=" + select_all(input_5)] else ["input_5=Empty"]
        File outputs_file = "outputs.txt"
    }
}