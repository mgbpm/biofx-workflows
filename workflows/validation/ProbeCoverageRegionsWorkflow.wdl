version 1.0

import "../../steps/DepthOfCoverage.wdl"
import "../../steps/CramToBamWorkflow.wdl"
import "../../steps/Utilities.wdl"

workflow ProbeCoverageRegions {
  input {
    File   coveragebed
    Array[RoiAndRefGeneFilePair]
           roigenes
    File   genenames
    File   bam_or_cram
    File?  bai
    File   referencefasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File   referenceindex = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File   referencedict  = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

    String gatk_image     = "broadinstitute/gatk3:3.7-0"
    String cov_image      = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/coverage:20230630"
  }

  String format = if      (sub(bam_or_cram, "^.*\\.[Bb][Aa][Mm]$",     "") == "") then "BAM"
                  else if (sub(bam_or_cram, "^.*\\.[Cc][Rr][Aa][Mm]$", "") == "") then "CRAM"
                  else "UNSUPPORTED"

  if (format == "UNSUPPORTED") {
    call Utilities.FailTask as UnsupportedFormatError {
      input:
        error_message = "the extension of ~{bam_or_cram} is neither .bam nor .cram"
    }
  }

  if (format == "BAM" && !defined(bai)) {
    call Utilities.FailTask as MissingBaiError {
      input:
        error_message = "bai is required when input is a BAM file"
    }
  }

  if (format == "CRAM" && defined(bai)) {
    call Utilities.FailTask as SpuriousBaiError {
      input:
        error_message = "bai must not be provided when input is a CRAM file"
    }
  }

  if (   (format == "BAM"  &&  defined(bai))
      || (format == "CRAM" && !defined(bai))) {

    if (format == "CRAM") {
      call CramToBamWorkflow.ConvertCramToBam {
        input:
            input_cram      = bam_or_cram
          , ref_fasta       = referencefasta
          , ref_fasta_index = referenceindex
          , ref_dict        = referencedict
      }
    }

    File effective_bam = select_first([ConvertCramToBam.output_bam, bam_or_cram])
    File effective_bai = select_first([ConvertCramToBam.output_bai, bai])

    call DepthOfCoverage.DepthOfCoverageWorkflow {
      input:
          run_wgs           = false
        , output_basename   = sub(basename(bam_or_cram), "\\.[^\\.]+$", "") + ".cov"
        , ref_fasta         = referencefasta
        , ref_fasta_index   = referenceindex
        , ref_dict          = referencedict
        , bam               = effective_bam
        , bai               = effective_bai
        , roi_all_bed       = coveragebed
        , roi_genes         = roigenes
        , gene_names        = genenames
        , cov_docker_image  = cov_image
        , gatk_docker_image = gatk_image
    }
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
