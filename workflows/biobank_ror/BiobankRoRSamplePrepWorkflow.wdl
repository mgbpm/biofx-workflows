version 1.0

import "../../steps/VCFUtils.wdl" # for processing input data
import "../../steps/Utilities.wdl" # for fail task
import "../../steps/FileUtils.wdl" # for fetching files

workflow IndividualSamplePrepWorkflow {
    input {
        # Dataset prep inputs
        Array[String] filenames
        Array[String] file_locations
        String fetch_file_type # either vcf or bcf
        File? sample_ids_list
        String dataset
        String dataset_structure # either "joint" or "individual"
        File? target_roi_bed
        Boolean is_sharded
        # Docker images
        String bcftools_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
        String orchutils_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230921"
        String ubuntu_docker_image = "ubuntu:latest"
        # GCP project and Terra workspace for secret retrieval
        String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
        String workspace_name
    }

    ## Test inputs and fail if not correct or compatible
    # Test the dataset_structure string parameter
    if ((dataset_structure != "individual") && (dataset_structure != "joint")) {
        call Utilities.FailTask as StructureFail {
            input:
                error_message = "The dataset_structure input parameter must be either 'individual' or 'joint'."
        }
    }
    # Test that the number of sample IDs matches the number of sample data locations
    if (length(filenames) != length(file_locations)) {
        call Utilities.FailTask as FileArrayLengthFail {
            input:
                error_message = "The number of dataset files does not match the number of file locations."
        }
    }
    # Test dataset structure input compatibility: joint files
    if (dataset_structure == "joint") {
        if (defined(sample_ids_list)){
            # Check if number of sample IDs is greater than 2
            if (length(read_lines(select_first([sample_ids_list]))) < 2) {
                call Utilities.FailTask as MultipleIDsFail {
                    input:
                        error_message = "There are not multiple IDs in the sample IDs list."
                }
            }
        }
        # Check if there are multiple files provided (if the dataset is defined as sharded)
        if ((is_sharded) && (length(filenames) < 2)) {
            call Utilities.FailTask as JointShardFail {
                input:
                    error_message = "Dataset structure is defined as joint and sharded, but multiple input files are not provided."
            }
        }
        # Check if there is only one file (if the dataset is not defined as sharded)
        if ((!is_sharded) && (length(filenames) > 1)) {
            call Utilities.FailTask as SingleJointVCFFail {
                input:
                    error_message = "Dataset structure is defined as joint and not sharded, but multiple input files are given."
            }
        }
    }
    # Test dataset structure input compatibility: individual samples
    if (dataset_structure == "individual") {
        if (defined(sample_ids_list)) {
            # Check if dataset has only one sample ID
            Array[String] individual_id = read_lines(select_first([sample_ids_list]))
            if (length(individual_id) > 1) {
                call Utilities.FailTask as IndividualIDFail {
                    input:
                        error_message = "There is more than one ID in the sample IDs list."
                }
            }
        }
        # Check if dataset is labeled as sharded; this check can be removed if datasets are actually sharded
        if (is_sharded) {
            call Utilities.FailTask as IndividualShardFail {
                input:
                    error_message = "Check dataset_structure and is_sharded input parameters. Individual sample VCFs/BCFs are not likely to be sharded."
            }
        }
    }

    ## Fetch the dataset's VCF/BCF files and their index files
    if ((fetch_file_type != "vcf") && (fetch_file_type != "bcf")) {
        call Utilities.FailTask as FetchFileTypeFail {
            input:
                error_message = "File type to fetch should be either 'vcf' or 'bcf'."
        }
    }
    scatter (i in range(length(filenames))) {
        call FileUtils.FetchFilesTask as FetchFiles {
            input:
                data_location = file_locations[i],
                recursive = true,
                file_types = [ fetch_file_type ],
                file_match_keys = [ filenames[i] ],
                docker_image = orchutils_docker_image,
                gcp_project_id = gcp_project_id,
                workspace_name = workspace_name
        }
        if (fetch_file_type == "vcf") {
            if (!defined(FetchFiles.vcf)) {
                call Utilities.FailTask as VCFNotFound {
                    input:
                        error_message = "VCF file " + filenames[i] + " not found in " + file_locations[i]
                }
            }
            if (!defined(FetchFiles.vcf_index)) {
                call Utilities.FailTask as VCFIndexNotFound {
                    input:
                        error_message = "Index for vcf " + filenames[i] + " not found in " + file_locations[i]
                }
            }
        }
        if (fetch_file_type == "bcf") {
            if (!defined(FetchFiles.bcf)) {
                call Utilities.FailTask as BCFNotFound {
                    input:
                        error_message = "BCF file " + filenames[i] + " not found in " + file_locations[i]
                }
            }
        }
    }

    ## Coerce Array[File?] to Array[File] for future tasks (and convert to VCF if necessary)
    if (fetch_file_type == "vcf") {
        Array[File] fetched_vcfs = select_all(FetchFiles.vcf)
        Array[File] fetched_vcfs_idx = select_all(FetchFiles.vcf_index)
    }
    if (fetch_file_type == "bcf") {
        Array[File] fetched_bcfs = select_all(FetchFiles.bcf)
    }
    Array[File] dataset_files = select_first([fetched_bcfs, fetched_vcfs])
    
    ## Create a sample IDs list (if one doesn't already exist)
    if (!defined(sample_ids_list)) {
        call FindSampleIDs {
            input:
                input_files = dataset_files,
                dataset = dataset,
                docker_image = bcftools_docker_image
        }
        # Check if the number of sample IDs matches dataset type input
        if ((length(read_lines(FindSampleIDs.sample_ids_list)) > 1) && dataset_structure == "individual") {
            call Utilities.FailTask as IndividualIDFail2 {
                input:
                    error_message = "There is more than one ID in the created sample IDs list."
            }
        }
        if ((length(read_lines(FindSampleIDs.sample_ids_list)) < 2) && dataset_structure == "joint" ) {
            call Utilities.FailTask as MultipleIDsFail2 {
                input:
                    error_message = "There are not multiple IDs in the created sample IDs list."
            }
        }
    }
    File sample_ids_file = select_first([sample_ids_list, FindSampleIDs.sample_ids_list])

    ## Prep data according to dataset structure inputs
    # Filter each VCF to regions in BED file if necessary
    if (defined(target_roi_bed)) {
        scatter (i in range(length(dataset_files))) {
            if (fetch_file_type == "vcf") {
                call VCFUtils.FilterVCFWithBEDTask as FilterVCFs {
                    input:
                        input_vcf = dataset_files[i],
                        input_vcf_idx = select_first([fetched_vcfs_idx])[i],
                        input_bed = select_first([target_roi_bed]),
                        regions_or_targets = "regions",
                        docker_image = bcftools_docker_image,
                        preemptible = 0
                }
                
            }
            if (fetch_file_type == "bcf") {
                call VCFUtils.FilterVCFWithBEDTask as FilterBCFs {
                    input:
                        input_vcf = dataset_files[i],
                        input_bed = select_first([target_roi_bed]),
                        regions_or_targets = "regions",
                        docker_image = bcftools_docker_image,
                        preemptible = 0
                }
            }
        }
        # Coerce filtered files Array[File?] to Array[File] for future tasks 
        if (fetch_file_type == "vcf") {
            Array[File] filtered_vcfs = select_all(FilterVCFs.output_vcf_gz)
        }
        if (fetch_file_type == "bcf") {
            Array[File] filtered_bcfs = select_all(FilterBCFs.output_vcf_gz)
        }
        Array[File] filtered_files = select_first([filtered_bcfs, filtered_vcfs])
    }
    # If data structure is joint, concat the shards (if necessary) and make individual sample VCFs
    if (dataset_structure == "joint") {
        if (is_sharded) {
            call VCFUtils.ConcatVCFsTask as ConcatJointVCFs {
                input:
                    input_vcfs = select_first([filtered_files, dataset_files]),
                    sorted = true,
                    output_basename = dataset + ".concat",
                    preemptible = 0,
                    docker_image = bcftools_docker_image
            }
        }
        call BatchSamplesTask as MakeBatches {
            input:
                sample_ids_list = sample_ids_file,
                docker_image = ubuntu_docker_image
        }
        scatter (samples_batch in MakeBatches.sample_batches) {
            Int individual_vcf_disk_space = ceil(length(read_lines(samples_batch)) / 10) + 50
            call MakeIndividualVCFsTask as MakeIndividualVCFs {
                input:
                    input_vcf = select_first([ConcatJointVCFs.output_vcf_gz, select_first([filtered_files, dataset_files])[0]]),
                    sample_ids_list = samples_batch,
                    output_basename = dataset,
                    docker_image = bcftools_docker_image,
                    disk_size = individual_vcf_disk_space
            }
            # Write a table containing individual VCFs info
            call WriteTSVTask as WriteIndividualVCFsTSV {
                input:
                    input_files = MakeIndividualVCFs.output_vcf_gz,
                    final_tsv = false,
                    dataset = dataset,
                    output_basename = dataset,
                    docker_image = bcftools_docker_image,
                    disk_size = individual_vcf_disk_space
            }
        }
        # Merge tables containing individual VCFs info (for uploading to Terra)
        call MergeTSVsTask as MergeIndividualVCFsTSVs {
            input:
                input_tsvs = WriteIndividualVCFsTSV.output_tsv,
                output_basename = dataset,
                docker_image = ubuntu_docker_image
        }
    }
    # If data structure is individual VCFs, concat them for Alamut annotation and loading
    if (dataset_structure == "individual") {
        call VCFUtils.ConcatVCFsTask as ConcatIndividualVCFs {
            input:
                input_vcfs = select_first([filtered_files, dataset_files]),
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
    if (!defined(MakeIndividualVCFs.output_vcf_gz)) {
        call WriteTSVTask as WriteOtherTSV {
            input:
                input_files = select_first([filtered_files, dataset_files]),
                final_tsv = true,
                dataset = dataset,
                output_basename = dataset,
                docker_image = bcftools_docker_image
        }
    }

    output {
        # Filtered VCF(s)/BCF(s)
        Array[File]? filtered_dataset_files = filtered_files
        # Concat VCF
        File concat_vcf = select_first([ConcatJointVCFs.output_vcf_gz, ConcatIndividualVCFs.output_vcf_gz])
        # Individual VCFs
        Array[File] individual_vcfs = flatten(select_first([MakeIndividualVCFs.output_vcf_gz]))
        # File to use for batch annotation
        File batch_annotation_input_file = select_first([ConcatJointVCFs.output_vcf_gz, ConcatIndividualVCFs.output_vcf_gz, select_first([FilterVCFs.output_vcf_gz, dataset_files])[0]])
        # List of paths for individual VCFs (Table to upload to Terra)
        File dataset_sample_table = select_first([MergeIndividualVCFsTSVs.output_tsv, WriteOtherTSV.output_tsv])
    }
}

task FindSampleIDs {
    input {
        Array[File] input_files
        String dataset
        String docker_image
        Int disk_size = ceil(size(input_files, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Gather IDs from all input files
        for file in '~{sep="' '" input_files}'
        do
            bcftools query -l $file >> "~{dataset}_ids_tmp.txt"
        done

        # Get only the unique IDs
        sort "~{dataset}_ids_tmp.txt" | uniq > "~{dataset}_ids.txt"
    >>>

   runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File sample_ids_list = "~{dataset}_ids.txt"
    }
}

task BatchSamplesTask {
    input {
        File sample_ids_list
        Int batch_size = 5000
        String docker_image
        Int disk_size = 10 + ceil(size(sample_ids_list, "GB"))
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        mkdir batches
        # Split the sample IDs list into batches of samples
        split "~{sample_ids_list}" "batches/batch_" -l "~{batch_size}" -d
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        Array[File] sample_batches = glob("batches/*")
    }
}

task MakeIndividualVCFsTask {
    input {
        File input_vcf
        File sample_ids_list
        String docker_image
        String output_basename
        Int disk_size
        Int mem_size = 4
        Int preemptible = 1
    }

     command <<<
        set -euxo pipefail

        # Create a file of IDs and output basename for bcftools split
            # Need 3 columns: ID, new ID, and output file base name
        while read -r line
        do
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
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        Array[File] output_vcf_gz = glob("output/*.vcf.gz")
    }
}

task WriteTSVTask {
    input {
        Array[File] input_files
        Boolean final_tsv
        String dataset
        String docker_image
        String output_basename
        Int disk_size = 10 + ceil(size(input_files, "GB"))
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # If the tsv is the final output tsv, then write the header first
            # Otherwise this will be written during the merging of all tsvs
        if [ "~{final_tsv}" == "true" ]
        then
            # dataset_sample_id: dataset and sample ID, e.g. Biobank 1004 subject 10000054 = “1004-10000054”
            # vcf_file: path to the individual VCF file for the sample
            # dataset_id: dataset number, e.g. Biobank 1004 = "1004"
            # sample_id: sample ID
            printf "entity:dataset_sample_id\tvcf_file\tdataset_id\tsample_id\n" > "~{output_basename}_dataset_sample_table.tsv"
        fi

        # Write info for each individual sample VCF to the dataset_sample_table.tsv
        for c in '~{sep="' '" input_files}'
        do
            vcf_path=$(echo $c | sed 's/\/cromwell_root\//gs:\/\//')
            sample_id=$(bcftools query --list-samples $c)
            dataset_sample_id="~{dataset}-${sample_id}"
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

task MergeTSVsTask {
    input {
        Array[File] input_tsvs
        String output_basename
        String docker_image
        Int disk_size = 10
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # Write the tsv header:
            # dataset_sample_id: dataset and sample ID, e.g. Biobank 1004 subject 10000054 = “1004-10000054”
            # vcf_file: path to the individual VCF file for the sample
            # dataset_id: dataset number, e.g. Biobank 1004 = "1004"
            # sample_id: sample ID
        printf "entity:dataset_sample_id\tvcf_file\tdataset_id\tsample_id\n" > "~{output_basename}_dataset_sample_table.tsv"

        # Add info from each tsv to the final output
        for c in '~{sep="' '" input_tsvs}'
        do
            cat $c >> "~{output_basename}_dataset_sample_table.tsv"
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