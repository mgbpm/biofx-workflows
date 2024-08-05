version 1.0

import "../../steps/FileUtils.wdl"
import "../../steps/Utilities.wdl"

## Requirements/expectations :
## - Two analysis-ready CRAM files for a single sample (as identified in RG:SM)
##
## Outputs :
## - One Merged CRAM file and its index


# WORKFLOW DEFINITION 
workflow MergeCrams {
  input {
    # GCP project and Terra workspace for secret retrieval
    String gcp_project_id = "mgb-lmm-gcp-infrast-1651079146"
    String workspace_name
    # Orchestration utils docker
    String orchutils_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
    # bcftools docker image
    String bcftools_docker_image = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/bcftools:1.17"
    Boolean fetch_files_verbose = false
    Int fetch_disk_size = 75
    Array[FileMatcher]? fetch_cram_file_matchers

    # Merge files Inputs
    String input_cram_1
    String input_cram_1_index
    String input_cram_2
    String input_cram_2_index
    String output_staging_bucket
    String output_sample_name

    # reference genome files
    String reference_build = "GRCh38"
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String basename_1 = basename(input_cram_1)
    String basename_2 = basename(input_cram_2)
    String sample_data_location_1 = sub(input_cram_1, "~{basename_1}", "")
    String sample_data_location_2 = sub(input_cram_2, "~{basename_2}", "")

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gatk_path = "/gatk/gatk"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String samtools_path = "samtools"
  }  

    #is the input a cram file?
    Boolean is_cram = sub(basename(input_cram_1), ".*\\.", "") == "cram"

    String sample_basename = if is_cram then  basename(input_cram_1, ".cram") else basename(input_cram_1, ".bam")
    String vcf_basename = sample_basename
    Array[String] fetch_cram_filter_keys = [sample_basename]

# fetch CRAM One
    call FileUtils.FetchFilesTask as FetchCramOne {
        input:
            data_location = sample_data_location_1,
            file_types = if defined(fetch_cram_file_matchers) then [] else ["cram"],
            recursive = false,
            file_match_keys = if defined(fetch_cram_file_matchers) then [] else fetch_cram_filter_keys,
            file_matchers = fetch_cram_file_matchers,
            verbose = fetch_files_verbose,
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            disk_size = fetch_disk_size
    }

    # Validate that we got workable files from the data location
    if (!defined(FetchCramOne.cram)) {
        call Utilities.FailTask as MissingBamOrCramFailure {
            input:
                error_message = "CRAM file not found in ~{sample_data_location_1}"
        }
    }
    if (defined(FetchCramOne.cram) && !defined(FetchCramOne.crai)) {
        call Utilities.FailTask as MissingCraiFailure {
            input:
                error_message = "Index file for CRAM " + basename(select_first([FetchCramOne.cram])) + " not found"
        }
    }

# fetch CRAM Two
    call FileUtils.FetchFilesTask as FetchCramTwo {
        input:
            data_location = sample_data_location_1,
            file_types = if defined(fetch_cram_file_matchers) then [] else ["cram"],
            recursive = false,
            file_match_keys = if defined(fetch_cram_file_matchers) then [] else fetch_cram_filter_keys,
            file_matchers = fetch_cram_file_matchers,
            verbose = fetch_files_verbose,
            docker_image = orchutils_docker_image,
            gcp_project_id = gcp_project_id,
            workspace_name = workspace_name,
            disk_size = fetch_disk_size
    }

    # Validate that we got workable files from the data location
    if (!defined(FetchCramTwo.cram)) {
        call Utilities.FailTask as MissingBamOrCramFailureTwo {
            input:
                error_message = "CRAM file not found in ~{sample_data_location_2}"
        }
    }
    if (defined(FetchCramTwo.cram) && !defined(FetchCramTwo.crai)) {
        call Utilities.FailTask as MissingCraiFailureTwo {
            input:
                error_message = "Index file for CRAM " + basename(select_first([FetchCramTwo.cram])) + " not found"
        }
    }

    call MergeCramsTask{
        input:
          input_cram_one = FetchCramOne.cram,
          input_cram_two = FetchCramTwo.cram,
          #sample_name = basename(input_cram_one),
          sample_name = sample_basename,
          output_sample_name = output_sample_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          docker = gitc_docker,
          samtools_path = samtools_path
    }

#push cram and index to gcp staging bucket
  #Transfer cram to staging bucket
  call FileUtils.CopyFilesTask as CopyCRAMToBucket {
  input:
      source_location = MergeCramsTask.output_cram,
      file_types = [".cram"],
      file_match_keys = [],
      file_matchers = [],
      target_location = output_staging_bucket,
      flatten = false,
      recursive = true,
      verbose = true,
      docker_image = orchutils_docker_image,
      disk_size = 75,
      gcp_project_id = gcp_project_id,
      workspace_name = workspace_name,
  }

  #Transfer gvcf index to staging bucket
  call FileUtils.CopyFilesTask as CopyCRAMIndexToBucket {
  input:
      source_location = MergeCramsTask.output_cram_index,
      file_types = [".cram.crai"],
      file_match_keys = [],
      file_matchers = [],
      target_location = output_staging_bucket,
      flatten = false,
      recursive = true,
      verbose = true,
      docker_image = orchutils_docker_image,
      disk_size = 75,
      gcp_project_id = gcp_project_id,
      workspace_name = workspace_name,
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_cram = MergeCramsTask.output_cram
    File output_cram_index = MergeCramsTask.output_cram_index
  }
}


# TASK DEFINITIONS

task MergeCramsTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? input_cram_one
    File? input_cram_two
    String sample_name
    String output_sample_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    String samtools_path
  }
    Float output_bam_one_size = size(input_cram_one, "GB") / 0.40
    Float output_bam_two_size = size(input_cram_two, "GB") / 0.40
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram_one, "GB") + size(input_cram_two, "GB") + output_bam_one_size + output_bam_two_size + ref_size + 20)
  
  command {
    set -e
    set -o pipefail

    echo "[`TZ="EST" date`] PWD: $PWD"

    #run cram one to bam conversion
    echo "[`TZ="EST" date`] converting cram 1 to bam"
    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram_one} |
    ~{samtools_path} view -b -o ~{sample_name}_one.bam
    ~{samtools_path} index -b ~{sample_name}_one.bam
    mv ~{sample_name}_one.bam.bai ~{sample_name}_one.bai
    echo "[`TZ="EST" date`] cram 1 to bam complete"

    #run cram two to bam conversion
    echo "[`TZ="EST" date`] converting cram 2 to bam"
    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram_two} |
    ~{samtools_path} view -b -o ~{sample_name}_two.bam
    ~{samtools_path} index -b ~{sample_name}_two.bam
    mv ~{sample_name}_two.bam.bai ~{sample_name}_two.bai
    echo "[`TZ="EST" date`] cram 2 to bam complete"

    #merge bams together
    echo "[`TZ="EST" date`] Merging bam 1: ~{sample_name}_one.bam"
    echo "[`TZ="EST" date`] With bam 2: ~{sample_name}_two.bam"
    echo "[`TZ="EST" date`] Writing merged output to: ~{output_sample_name}.bam"
    ~{samtools_path} merge --threads 4 -f ~{output_sample_name}.bam ~{sample_name}_one.bam ~{sample_name}_two.bam 
    echo "[`TZ="EST" date`] bam merge complete"

    #convert merged BAMS back to CRAM and generate index
    echo "[`TZ="EST" date`] converting merged bam to cram"
    ~{samtools_path} view -C -T ~{ref_fasta} -o ~{output_sample_name}.cram ~{output_sample_name}.bam
    ~{samtools_path} index -@ 4 -c ~{output_sample_name}.cram
    #mv ~{output_sample_name}.cram.crai ~{output_sample_name}.crai
    echo "[`TZ="EST" date`] ran merged bam to cram"

  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 2])
 }
  output {
    File output_cram = "~{output_sample_name}.cram"
    File output_cram_index = "~{output_sample_name}.crai"
  }
}
