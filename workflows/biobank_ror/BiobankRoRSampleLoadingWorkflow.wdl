version 1.0

import "../../steps/QCEval.wdl" # for qc of each vcf
import "../../steps/FASTUtils.wdl" # for loading to FAST

workflow SampleLoadingWorkflow {
    input {
        # Dataset and sample inputs
        File input_vcf
        String sample_ID
        String dataset
        Boolean remove_extra_ad_values = false
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
        # Ubuntu docker image
        String ubuntu_docker_image = "ubuntu:latest"
        # Orchestration utils docker image
        String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20231129"
        # reference genome
        String reference_build = "GRCh38"
        # qceval inputs
        String qceval_project_type
        String qceval_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/qceval:20231005"
        # FAST loading inputs
        Boolean has_haploid_sites = false
        String sample_data_load_config_name = "Sample_VCF_PPM_Eval"
    }

    # Remove extra values in AD field if necessary
    if (remove_extra_ad_values) {
        call RemoveExtraADValuesTask as RemoveExtraAD {
            input:
                input_vcf = input_vcf,
                docker_image = ubuntu_docker_image
        }
    }
    # Annotate individual sample VCF with QCEval INFO field
    call QCEval.QCEvalTask {
        input:
            input_vcf = select_first([RemoveExtraAD.fixed_vcf, input_vcf]),
            project_type = qceval_project_type,
            output_basename = dataset + "_" + sample_ID + ".qceval",
            docker_image = qceval_docker_image
    }
    # Load sample data to FAST
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

    output {
        # "Fixed" AD field VCF
        File? fixed_ad_field_vcf = RemoveExtraAD.fixed_vcf
        # List of variants with removed AD values
        File? fixed_variants_list = RemoveExtraAD.fixed_variants_list
        # Annotated qceval VCF
        File qceval_vcf_gz = QCEvalTask.output_vcf_gz
        # FAST sample data name
        String fast_sample_data_name = dataset + "_" + sample_ID
        # QCEval load data id
        String qceval_data_load_id = QCEvalLoadTask.data_load_id
    }
}

task RemoveExtraADValuesTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String docker_image
        Int disk_size = ceil(size(input_vcf, "GB") * 2) + 5
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        mkdir --parents OUTPUT
        mkdir --parents WORK

        # Prep VCF
        cp "~{input_vcf}" "WORK/~{output_basename}.vcf.gz"
        gunzip "WORK/~{output_basename}.vcf.gz"
        grep "#" "WORK/~{output_basename}.vcf" > "OUTPUT/~{output_basename}.fixed_values.vcf"
        grep -v "#" "WORK/~{output_basename}.vcf" > "WORK/~{output_basename}.no_header.vcf" 
    
        while read -r line; do
            # Find the values in the AD field
            ad_values=$( echo "$line" | awk '{print $10}' | cut -d ":" -f 2 )
            tmp_num_ad_values=$( echo "$ad_values" | grep -o "," | wc -l )
            num_ad_values=$( expr $tmp_num_ad_values + 1 ) 

            # Find the number of reference alleles
            tmp_num_ref_values=$( echo "$line" | awk '{print $4}' | grep -o "," | wc -l || : )
            num_ref_values=$( expr $tmp_num_ref_values + 1 )
            
            # Find the number of alternate alleles
            tmp_num_alt_values=$( echo "$line" | awk '{print $5}' | grep -o "," | wc -l || : )
            num_alt_values=$( expr $tmp_num_alt_values + 1 )
            
            # Find the number of expected values in the AD field
            expected_num_ad_values=$( expr $num_ref_values + $num_alt_values )
    
            # If the number of AD values is greater than expected
                # Otherwise, print the unchanged line to the new VCF
            if [[ $num_ad_values -gt $expected_num_ad_values ]]; then
                # Log the variant as being "fixed"
                echo "$line" | awk '{print $1 ":" $2 "_" $4 "_" $5}' >> "OUTPUT/~{output_basename}_fixed_variants.txt"

                # Find if the last value in the AD field is a zero
                if [[ $( echo "$ad_values" | cut -d "," -f $num_ad_values ) -ne 0 ]]; then
                    echo "ERROR: Last value in AD field was not a zero" 1>&2
                    exit 1
                else
                    new_ad_vals=$( echo "$ad_values" | cut -d "," -f $num_ad_values --complement )
                    new_line=$( echo "$line" | sed 's/'$ad_values'/'$new_ad_vals'/' )
                    printf "$new_line\n" >> "OUTPUT/~{output_basename}.fixed_values.vcf"
                fi
            else
                echo "$line" >> "OUTPUT/~{output_basename}.fixed_values.vcf"
            fi
        done < "WORK/~{output_basename}.no_header.vcf"

        # Zip the new VCF
        gzip "OUTPUT/~{output_basename}.fixed_values.vcf"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File fixed_vcf = "OUTPUT/~{output_basename}.fixed_values.vcf"
        File fixed_variants_list = "OUTPUT/~{output_basename}_fixed_variants.txt"
    }
}