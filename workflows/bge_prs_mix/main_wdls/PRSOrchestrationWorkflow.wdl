version 1.0

import "../../../steps/Utilities.wdl"
import "../../lowpassimputation/Glimpse2Imputation.wdl"
import "../subwdls/PRSRawScoreWorkflow.wdl"
import "../subwdls/PRSMixScoreWorkflow.wdl"
import "../subwdls/PRSPCAWorkflow.wdl"
import "../tasks/PRSStructs.wdl"

workflow PRSOrchestrationWorkflow {
    input {
        # GLIMPSE2 INPUTS
        File? glimpse_reference_chunks
        Array[File]? input_crams
        Array[File]? input_crai
        Array[String]? sample_ids
        String? vcf_basename
        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        Boolean collect_glimpse_qc = true
        String glimpse_docker_image = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_e0b9b56"
        String glimpse_extract_docker_image = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        File? glimpse_monitoring_script
        File? ref_fasta
        File? ref_fai
        File? ref_dict

        # PRS INPUTS
        Array[File] condition_model_manifests
        File condition_yaml
        
        # DEBUGGING INPUTS
        Boolean run_glimpse = true
        File? query_vcf
    }

    # Run input checks
    if (run_glimpse) {
        if (!defined(glimpse_reference_chunks) || !defined(input_crams) || !defined(input_crai)) {
            String GlimpseInputError = "Missing one or more GLIMPSE inputs: ref chunks, input crams, and/or input crai."
        }
        if (!defined(ref_fasta) || !defined(ref_fai) || !defined(ref_dict)) {
            String GlimpseReferenceInputError = "Missing one or more reference files for GLIMPSE: fasta, fai, and/or dict."
        }
    }

    if (!run_glimpse && !defined(query_vcf)) {
        String InputVcfError = "If GLIMPSE is not being run, please input a VCF for running PRS modules."
    }

    String input_check_result = select_first([GlimpseInputError, GlimpseReferenceInputError, InputVcfError, "No error"])

    if (input_check_result != "No error") {
        call Utilities.FailTask as FailInputCheck {
            input:
                error_message = input_check_result
        }
    }

    if (input_check_result == "No error") {
        if (run_glimpse) {
            # Run GLIMPSE to get imputed low-pass variants
            call Glimpse2Imputation.Glimpse2Imputation as RunGlimpse {
                input:
                    reference_chunks = select_first([glimpse_reference_chunks]),
                    crams = select_first([input_crams]),
                    cram_indices = select_first([input_crai]),
                    sample_ids = select_first([sample_ids]),
                    fasta = select_first([ref_fasta]),
                    fasta_index = select_first([ref_fai]),
                    ref_dict = select_first([ref_dict]),
                    output_basename = select_first([vcf_basename, "PRS_Glimpse"]),
                    impute_reference_only_variants = impute_reference_only_variants,
                    call_indels = call_indels,
                    n_burnin = select_first([n_burnin]),
                    n_main = select_first([n_main]),
                    effective_population_size = select_first([effective_population_size]),
                    collect_qc_metrics = collect_glimpse_qc,
                    preemptible = 9,
                    docker = glimpse_docker_image,
                    docker_extract_num_sites_from_reference_chunk = glimpse_extract_docker_image,
                    cpu_ligate = 4,
                    mem_gb_ligate = 4,
                    monitoring_script = select_first([glimpse_monitoring_script])
            }
        }

        scatter (i in range(length(condition_model_manifests))) {
            String condition_name = sub(basename(condition_model_manifests[i]), "_[[0-9]]+\\.json$", "")

            AdjustmentModelData model_data = read_json(condition_model_manifests[i])

            # Get PRS raw scores for each condition
            call PRSRawScoreWorkflow.PRSRawScoreWorkflow as PRSRawScores {
                input:
                    input_vcf = select_first([RunGlimpse.imputed_vcf, query_vcf]),
                    adjustment_model_manifest = condition_model_manifests[i]
            }

            # Get the PRS mix raw score for each condition
            call PRSMixScoreWorkflow.PRSMixScoreWorkflow as PRSMixScores {
                input:
                    condition_name = condition_name,
                    raw_scores = PRSRawScores.prs_raw_scores,
                    score_weights = select_first([model_data.score_weights])
            }

            # Perform PCA with population model
            call PRSPCAWorkflow.PRSPCAWorkflow as PerformPCA {
                input:
                    condition_name = condition_name,
                    input_vcf = select_first([RunGlimpse.imputed_vcf, query_vcf]),
                    adjustment_model_manifest = condition_model_manifests[i],
                    prs_raw_scores = PRSMixScores.prs_mix_raw_score
            }
        }

        # Categorize each condition's score into bins; report percentile & bin
        call SummarizeScores {
            input:
                scores = select_first([select_all(PerformPCA.adjusted_scores)]),
                condition_yaml = condition_yaml
        }
    }

    output {
        # Glimpse outputs
        File? glimpse_vcf = RunGlimpse.imputed_vcf
        File? glimpse_vcf_index = RunGlimpse.imputed_vcf_index
        File? glimpse_qc_metrics = RunGlimpse.qc_metrics
        Array[File?]? glimpse_phase_monitoring = RunGlimpse.glimpse_phase_monitoring
        File? glimpse_ligate_monitoring = RunGlimpse.glimpse_ligate_monitoring

        # PRS Outputs
        Array[Array[File]]? prs_raw_scores = PRSRawScores.prs_raw_scores
        Array[File]? prs_mix_raw_score = PRSMixScores.prs_mix_raw_score
        Array[File?]? prs_adjusted_score = PerformPCA.adjusted_scores
        Array[File]? pc_projection = PerformPCA.pc_projection
        Array[File]? pc_plot = PerformPCA.pc_plot

        # Individual Outputs
        Array[File]? individuals_risk_summaries = SummarizeScores.individual_risk_summaries
    }
}

task SummarizeScores {
    input {
        Array[File] scores
        File condition_yaml
        String docker_image = "python:3.11"
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