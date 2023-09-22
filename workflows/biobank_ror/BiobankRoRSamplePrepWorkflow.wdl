version 1.0

import "../../steps/VCFUtils.wdl" # for processing input data
import "../../steps/Utilities.wdl" # for fail task

workflow IndividualSamplePrepWorkflow {
    input {
        # Dataset prep inputs
        Array[Pair[File,File]] dataset_files
        File sample_ids_list
        String dataset
        String dataset_structure # either "joint" or "individual"
        File? target_roi_bed
        Boolean is_sharded
        # bcftools docker image
        String bcftools_docker_image
    }

    ## Test inputs and fail if not correct or compatible
    # Test the dataset_structure string parameter
    if ((dataset_structure != "individual") && (dataset_structure != "joint")) {
        call Utilities.FailTask as StructureFail {
            input:
                error_message = "The dataset_structure input parameter must be either 'individual' or 'joint'."
        }
    }
    # Test dataset structure input compatibility: joint VCFs
    if (dataset_structure == "joint") {
        Array[String] joint_ids = read_lines(sample_ids_list)
        # Check if number of sample IDs is greater than 2
        if (length(joint_ids) < 2) {
            call Utilities.FailTask as MultipleIDsFail {
                input:
                    error_message = "There are not multiple IDs in the sample IDs list."
            }
        }
        # Check if there are multiple shards provided (if the dataset is defined as sharded)
        if ((is_sharded) && (length(dataset_files) < 2)) {
            call Utilities.FailTask as JointShardFail {
                input:
                    error_message = "Dataset structure is defined as joint VCFs and sharded, but multiple VCFs are not provided in dataset_files input parameter."
            }
        }
        # Check if there is only one VCF (if the dataset is not defined as sharded)
        if ((!is_sharded) && (length(dataset_files) > 1)) {
            call Utilities.FailTask as SingleJointVCFFail {
                input:
                    error_message = "Dataset structure is defined as joint VCFs and not sharded, but multiple VCFs are given."
            }
        }
    }
    # Test dataset structure input compatibility: individual sample VCFs
    if (dataset_structure == "individual") {
        Array[String] individual_id = read_lines(sample_ids_list)
        # Check if dataset has only one sample ID
        if (length(individual_id) > 1) {
            call Utilities.FailTask as IndividualIDFail {
                input:
                    error_message = "There is more than one ID in the sample IDs list."
            }
        }
        # Check if dataset is labeled as sharded; this check can be removed if datasets are actually sharded
        if (is_sharded) {
            call Utilities.FailTask as IndividualShardFail {
                input:
                    error_message = "Check dataset_structure and is_sharded input parameters. Individual sample vcfs are not likely to be sharded."
            }
        }
    }
    # Test that the pairs of files in dataset files array match (e.g. that the VCF and index files in each pair are the same)
    scatter (pair in dataset_files) {
        File dataset_vcfs = pair.left
        File dataset_tbis = pair.right
        String basename = sub(basename(pair.left), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String basename_vcf = basename(dataset_vcfs)
        String basename_tbi = basename(dataset_tbis)
        # Check that VCFs are on the left and index files are on the right
        if (((basename + ".vcf.gz.tbi") == basename_vcf) || ((basename + ".vcf.gz") == basename_tbi)) {
            call Utilities.FailTask as ArrayConstructionFail {
                input:
                    error_message = "Ensure dataset VCF files are the left of the array pairs. Dataset VCF index files are the right of the array pairs."
            }
        }
        # Check that the index file is for the VCF in the pair
        if ((basename_vcf + ".tbi") != basename_tbi) {
            call Utilities.FailTask as MatchFail {
                input:
                    error_message = "A pair of input tbi and vcf names do not match. Check the pairs in dataset_files input parameter."
            }
        }
    }

    ## Prep data according to dataset structure inputs
    # Filter each VCF to regions in BED file if necessary
    if (defined(target_roi_bed)) {
        scatter (pair in dataset_files) {
            call VCFUtils.FilterVCFWithBEDTask as FilterVCFs {
                input:
                    input_vcf = pair.left,
                    input_vcf_idx = pair.right,
                    input_bed = select_first([target_roi_bed]),
                    regions_or_targets = "regions",
                    docker_image = bcftools_docker_image,
                    preemptible = 0
            }
        }
    }
    # If data structure is joint VCF, concat the shards (if necessary) and make individual VCFs
    if (dataset_structure == "joint") {
        if (is_sharded) {
            call VCFUtils.ConcatVCFsTask as ConcatJointVCFs {
                input:
                    input_vcfs = select_first([FilterVCFs.output_vcf_gz, dataset_vcfs]),
                    sorted = true,
                    output_basename = dataset + "_" + ".concat",
                    preemptible = 0,
                    docker_image = bcftools_docker_image
            }
        }
        call MakeIndividualVCFsTask as MakeIndividualVCFs {
            input:
                input_vcf = select_first([ConcatJointVCFs.output_vcf_gz, select_first([FilterVCFs.output_vcf_gz, dataset_vcfs])[0]]),
                sample_ids_list = sample_ids_list,
                output_basename = dataset,
                docker_image = bcftools_docker_image
        }
    }
    # If data structure is individual VCFs, concat them for Alamut annotation and loading
    if (dataset_structure == "individual") {
        call VCFUtils.ConcatVCFsTask as ConcatIndividualVCFs {
            input:
                input_vcfs = select_first([FilterVCFs.output_vcf_gz, dataset_vcfs]),
                sorted = true,
                output_basename = dataset + ".concat",
                preemptible = 0,
                docker_image = bcftools_docker_image
        }
        # If dataset is already individual VCFs, they are unlikely to be sharded. But if this occurs, the
        # solution can be developed here:
        # if ((dataset_structure == "individual") && (is_sharded)) {...}
        # NOTE: should remove the failure in the test section above if this is developed
    }

    ## Write a table containing individual VCFs info to upload to the Terra workspace
    Array[File] samples_to_load = select_first([MakeIndividualVCFs.output_vcf_gz, FilterVCFs.output_vcf_gz, dataset_vcfs])
    call WriteTsvTask as WriteTsv {
        input:
            input_files = samples_to_load,
            dataset = dataset,
            output_basename = dataset,
            docker_image = bcftools_docker_image
    }

    output {
        # Filtered VCF(s)
        Array[File]? filtered_vcfs = FilterVCFs.output_vcf_gz
        # Concat VCF
        File concat_vcf = select_first([ConcatJointVCFs.output_vcf_gz, ConcatIndividualVCFs.output_vcf_gz])
        # Individual VCFs
        Array[File]? individual_vcfs = MakeIndividualVCFs.output_vcf_gz
        # File to use for batch annotation
        File batch_annotation_input_file = select_first([ConcatJointVCFs.output_vcf_gz, ConcatIndividualVCFs.output_vcf_gz, select_first([FilterVCFs.output_vcf_gz, dataset_vcfs])[0]])
        # List of paths for individual VCFs (Table to upload to Terra)
        File dataset_sample_table = WriteTsv.output_tsv
    }
}

task MakeIndividualVCFsTask {
    input {
        File input_vcf
        File sample_ids_list
        String docker_image
        String output_basename
        Int disk_size = 10 + ceil(size(input_vcf, "GB") * 2)
        Int preemptible = 1
    }

     command <<<
        set -euxo pipefail

        # Create a file of IDs and output basename for bcftools split
        while read -r line
        do
            # Need 3 columns: ID, new ID, and output file base name
            printf "${line}\t${line}\t~{output_basename}_${line}\n" >> "~{output_basename}_sample_ids.txt"
        done < "~{sample_ids_list}"

        # Put all individual sample vcfs in an output directory
        mkdir "output"
        # Split the input vcf by sample
        bcftools +split -S "~{output_basename}_sample_ids.txt" --output-type z --output "output" "~{input_vcf}"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: preemptible
    }

    output {
        Array[File] output_vcf_gz = glob("output/*.vcf.gz")
    }
}

task WriteTsvTask {
    input {
        Array[File] input_files
        String dataset
        String output_basename
        String docker_image
        Int disk_size = 10 + ceil(size(input_files[0], "GB") * 2.2)
        Int preemptible = 1
    }

    command <<<
        # Write the tsv header:
            # dataset_sample_id: dataset and sample ID, e.g. Biobank 1004 subject 10000054 = “1004-10000054”
            # vcf_file: path to the individual VCF file for the sample
            # dataset_id: dataset number, e.g. Biobank 1004 = "1004"
            # sample_id: sample ID
        printf "entity:dataset_sample_id\tvcf_file\tdataset_id\tsample_id\n" > "~{output_basename}_dataset_sample_table.tsv"

        # Write info for each individual sample VCF to the dataset_sample_table.tsv
        for c in '~{sep="' '" input_files}'
        do
            # Find bucket location of VCF
            vcf_path=$(echo $c | sed 's/\/cromwell_root\//gs:\/\//')
            # Find the sample ID in the VCF
            sample_id=$(bcftools query --list-samples $c)
            # Find the dataset-sample ID
            dataset_sample_id="~{dataset}-${sample_id}"
            # Write the info
            printf "${dataset_sample_id}\t${vcf_path}\t~{dataset}\t${sample_id}\n" >> "~{output_basename}_dataset_sample_table.tsv"
        done
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        File output_tsv = "~{output_basename}_dataset_sample_table.tsv"
    }
}