version 1.0

## Copyright Broad Institute, 2019
## 
## The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
## from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
## When executed the workflow scatters the HaplotypeCaller tool over a sample
## using an intervals list file. The output file produced will be a
## single gvcf file which can be used by the joint-discovery workflow.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support 
## - Successfully tested on v53
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION 
workflow HaplotypeCallerGvcf_GATK4 {
  input {
    File input_bam
    File input_bam_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File scattered_calling_intervals_list
  
    Boolean make_gvcf = true
    Boolean make_bamout = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gatk_path = "/gatk/gatk"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String samtools_path = "samtools"
  }  

    Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

    #is the input a cram file?
    Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"

    String sample_basename = if is_cram then  basename(input_bam, ".cram") else basename(input_bam, ".bam")
    String vcf_basename = sample_basename
    String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String output_filename = vcf_basename + output_suffix

    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Int potential_hc_divisor = length(scattered_calling_intervals) - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  if ( is_cram ) {
    call CramToBamTask {
      input:
        input_cram = input_bam,
        sample_name = sample_basename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gitc_docker,
        samtools_path = samtools_path
    }
  }

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = select_first([CramToBamTask.output_bam, input_bam]),
        input_bam_index = select_first([CramToBamTask.output_bai, input_bam_index]),
        interval_list = interval_file,
        output_filename = output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_gvcf = make_gvcf,
        make_bamout = make_bamout,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_filename = output_filename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MergeGVCFs.output_vcf
    File output_vcf_index = MergeGVCFs.output_vcf_index
  }
}

# TASK DEFINITIONS

task CramToBamTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    String samtools_path
  }
    Float output_bam_size = size(input_cram, "GB") / 0.40
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20
  
  command {
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 2])
 }
  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    File interval_list
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    Int hc_scatter

    String? gcs_project_for_requester_pays

    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
  }

  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"]) 

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  String vcf_basename = if make_gvcf then  basename(output_filename, ".gvcf") else basename(output_filename, ".vcf")
  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      description: "a bam file",
      localization_optional: true
    }
    input_bam_index: {
      description: "an index file for the bam input",
      localization_optional: true
    }
  }
    command {
    set -e
  
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_filename} \
      -contamination ~{default="0" contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      --disable-spanning-event-genotyping true \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 2])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
    File bamout = "~{vcf_basename}.bamout.bam"
  }
}
# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1
  
  command {
  set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 2])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}

task GenomePanelsVariantCallingTask {
  input {
        # Command parameters
        File input_bam
        File input_bam_index
        File interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File dbsnp
        File dbsnp_vcf_index
        String gatk_path = "/gatk/gatk"

        # Runtime parameters
        String docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        Int? mem_gb
        Int? disk_space_gb
        Boolean use_ssd = false
        Int? preemptible_attempts
    }

    Int machine_mem_gb = select_first([mem_gb, 24])
    Int command_mem_gb = machine_mem_gb - 1

    Float ref_size = size(ref_fasta, "GB") + size(ref_dict, "GB") + size(dbsnp, "GB")
    Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

    Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"
    String sample_basename = if is_cram then  basename(input_bam, ".cram") else basename(input_bam, ".bam")
    String gvcf = sample_basename + ".g.vcf"
    String all_calls_vcf = sample_basename + ".allcalls.vcf"

    command <<<
      set -euxo pipefail

      ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G" \
      HaplotypeCaller \
      --input ~{input_bam} \
      --output ~{gvcf} \
      --reference ~{ref_fasta} \
      --intervals ~{interval_list} \
      --dbsnp ~{dbsnp} \
      -mbq 10 \
      -stand-call-conf 30 \
      -A QualByDepth \
      -A FisherStrand \
      -A RMSMappingQuality \
      -A MappingQualityZero \
      -A StrandOddsRatio \
      -A DepthPerAlleleBySample \
      --read-filter MappingQualityReadFilter \
      --minimum-mapping-quality '17' \
      --read-filter MappingQualityNotZeroReadFilter \
      -ip '15' \
      --do-not-run-physical-phasing 'true' \
      -ERC BP_RESOLUTION
      
      ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G" \
      GenotypeGVCFs \
      --variant ~{gvcf} \
      --output ~{all_calls_vcf} \
      --reference ~{ref_fasta} \
      --intervals ~{interval_list} \
      --dbsnp ~{dbsnp} \
      -A QualByDepth \
      -A FisherStrand \
      -A RMSMappingQuality \
      -A MappingQualityZero \
      -A StrandOddsRatio \
      -A DepthPerAlleleBySample \
      --read-filter MappingQualityReadFilter \
      --minimum-mapping-quality '17' \
      -ip '15' \
      --read-filter MappingQualityNotZeroReadFilter 

    >>>

    runtime {
      docker: docker
      memory: machine_mem_gb + " GB"
      disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
      preemptible: select_first([preemptible_attempts, 2])
    }

    output {
      File all_calls_vcf_file = "~{all_calls_vcf}"
      File all_calls_vcf_idx_file = "~{all_calls_vcf}.idx"
      File gvcf_file = "~{gvcf}"
      File gvcf_idx_file = "~{gvcf}.idx"
    }
}

task GenomePanelsRefSitesSortTask {
  input {
        # Command parameters
        File gvcf_file
        File all_calls_vcf_file
        File ref_positions_vcf_file
        String gatk_path = "/gatk/gatk"
        String java_path = "/usr/lib/jvm/java-8-openjdk-amd64/bin/java"

        # Runtime parameters
        String docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        Int? mem_gb
        Int? disk_space_gb
        Boolean use_ssd = false
        Int? preemptible_attempts
    }

    Int machine_mem_gb = select_first([mem_gb, 24])
    Int command_mem_gb = machine_mem_gb - 1

    Int disk_size = ceil(size(gvcf_file, "GB") + size(all_calls_vcf_file, "GB")) + 10

    String sample_basename = basename(gvcf_file, ".g.vcf")
    String all_bases_vcf = sample_basename + ".allbases.vcf"

    command <<<
      set -euxo pipefail

      ~{java_path} -Xms12g "-Xmx~{command_mem_gb}G" \
      -jar ~{gatk_path}.jar SortVcf \
      -I ~{ref_positions_vcf_file} \
      -I ~{all_calls_vcf_file} \
      -O ~{all_bases_vcf}

    >>>

    runtime {
      docker: docker
      memory: machine_mem_gb + " GB"
      disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
      preemptible: select_first([preemptible_attempts, 2])
    }

    output {
      File all_bases_vcf_file = "~{all_bases_vcf}"
      File all_bases_vcf_idx_file = "~{all_bases_vcf}.idx"
    }
}
