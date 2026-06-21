version 1.0

import "../../steps/DepthOfCoverage.wdl"

workflow ProbeCoverageRegions {
  input {
    File   coveragebed
    Array[RoiAndRefGeneFilePair]
           roigenes
    File   genenames      = "gs://lmm-reference-data/roi/HGNC_genenames_05272022.txt"
    File   inputcram
    File   referencefasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File   referenceindex = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File   referencedict  = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

    String samtools_image = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String gatk_image     = "broadinstitute/gatk3:3.7-0"
    String cov_image      = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/coverage:20230630"
  }

  call CramToBam {
    input:
        input_cram      = inputcram
      , sample_name     = basename(inputcram, ".cram")
      , ref_fasta       = referencefasta
      , ref_fasta_index = referenceindex
      , ref_dict        = referencedict
      , docker          = samtools_image
  }

  call DepthOfCoverage.DepthOfCoverageWorkflow {
    input:
        run_wgs           = false
      , output_basename   = basename(inputcram, ".cram") + ".cov"
      , ref_fasta         = referencefasta
      , ref_fasta_index   = referenceindex
      , ref_dict          = referencedict
      , bam               = CramToBam.output_bam
      , bai               = CramToBam.output_bai
      , roi_all_bed       = coveragebed
      , roi_genes         = roigenes
      , gene_names        = genenames
      , cov_docker_image  = cov_image
      , gatk_docker_image = gatk_image
  }

  output {                                                                                                                                # FROM:
    File  roi_sample_interval_summary                = select_first([DepthOfCoverageWorkflow.roi_sample_interval_summary])                # DepthOfCoverageROITask
    File  roi_sample_interval_statistics             = select_first([DepthOfCoverageWorkflow.roi_sample_interval_statistics])             # DepthOfCoverageROITask
    File  roi_sample_statistics                      = select_first([DepthOfCoverageWorkflow.roi_sample_statistics])                      # DepthOfCoverageROITask
    File  roi_sample_summary                         = select_first([DepthOfCoverageWorkflow.roi_sample_summary])                         # DepthOfCoverageROITask
    File  roi_sample_cumulative_coverage_counts      = select_first([DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_counts])      # DepthOfCoverageROITask
    File  roi_sample_cumulative_coverage_proportions = select_first([DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_proportions]) # DepthOfCoverageROITask

    File? mt_summary                                 = DepthOfCoverageWorkflow.mt_summary                                                 # DepthOfCoverageSummaryTask
    File  gene_summary                               = select_first([DepthOfCoverageWorkflow.gene_summary])                               # DepthOfCoverageSummaryTask
    File  gene_summary_unknown                       = select_first([DepthOfCoverageWorkflow.gene_summary_unknown])                       # DepthOfCoverageSummaryTask
    File  gene_summary_entrez                        = select_first([DepthOfCoverageWorkflow.gene_summary_entrez])                        # DepthOfCoverageSummaryTask
  }
}

task CramToBam {
  input {
    File   ref_fasta
    File   ref_fasta_index
    File   ref_dict
    File   input_cram
    String sample_name
    String docker
    Int    preemptible = 2
  }

  Float output_bam_size = size(input_cram, "GB") / 0.40
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int   disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20

  command <<<
  set -o errexit
  set -o pipefail

  samtools view -h -T '~{ref_fasta}' '~{input_cram}' |
  samtools view -b -o '~{sample_name}.bam' -
  samtools index -b '~{sample_name}.bam'
  mv '~{sample_name}.bam.bai' '~{sample_name}.bai'
  >>>

  runtime {
    docker      : docker
    memory      : "15 GB"
    disks       : "local-disk ~{disk_size} HDD"
    preemptible : preemptible
  }

  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}
