version 1.0

workflow PRSSummaryWorkflow {
    input {
        # PRS inputs
        Array[File] prs_scores
        File condition_yaml
        # Ubuntu Docker image
        String python_docker_image = "python:3.11"
    }

    # Get summary for each individual of all their PRS scores (one from each condition)
    call SummarizeScores {
        input:
            scores = prs_scores,
            condition_yaml = condition_yaml,
            docker_image = python_docker_image
    }

    output {
        Array[File] individual_risk_summaries = SummarizeScores.individual_risk_summaries
    }
}

task SummarizeScores {
    input {
        Array[File] scores
        File condition_yaml
        String docker_image
        Int disk_size = ceil(size(scores, "GB") + size(condition_yaml, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        pip install PyYAML

        mkdir -p WORK/summaries
        mkdir -p OUTPUT

        # Extract all sample IDs from a score file
            # NOTE: IDs should be the same across score files
        score_file_array=('~{sep="' '" scores}')
        sed '1d;' ${score_file_array[0]} | awk '{ print $2 }' > WORK/sample_ids.txt

        # For each scores file, summarize the file with the sample info, condition name, and bin count
        for c in '~{sep="' '" scores}'; do
            file_basename=$(basename $c .tsv)

            python3 -c '
import yaml

# Load configs for conditions/diseases
with open("~{condition_yaml}", "r") as yml_file:
    conditions_configs = yaml.safe_load(yml_file)

with open("'$c'", "r") as scores_file:
    condition_summary_file = open("WORK/summaries/'$file_basename'.summary.csv", "w")
    condition_summary_file.write("Sample,Condition,Risk,Bin_Count\n")

    header = scores_file.readline().strip("\n").replace(" ", "").split("\t")
    line = scores_file.readline().strip("\n").replace(" ", "").split("\t")

    while line != [""]:
        sample_id = line[header.index("IID")]

        percentile = str(line[header.index("percentile")])
        bins = str(conditions_configs["'$file_basename'"]["bin_count"])
        print(sample_id, percentile, bins)

        condition_summary_file.write(
                ",".join([
                    sample_id,
                    "'$file_basename'",
                    percentile,
                    bins
                ]) + "\n"
        )
            
        line = scores_file.readline().strip("\n").replace(" ", "").split("\t")

    condition_summary_file.close()'
        done
    
        # Get per sample summaries
        while read line; do
            printf "Condition,Risk,Bin_Count\n" >> OUTPUT/"${line}"_prs_summary.csv
            for file in WORK/summaries/*; do
                grep "${line}" $file | cut -d "," -f 2,3,4 >> OUTPUT/"${line}"_prs_summary.csv
            done
        done < WORK/sample_ids.txt
        
    >>>

    runtime {
          docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        Array[File] individual_risk_summaries = glob("OUTPUT/*")
    }
}