version 1.0
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# -----------------------------------------------------------------------------

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/steps/DepthOfCoverage.wdl"

workflow BgeMonogenicValidation {
  input {
    File   inputcram
    File   coveragebed
    Array[RoiAndRefGeneFilePair] roigenes
    File   referencefasta
    File   referenceindex
    File   referencedict
    File?  final_referencefasta
    File?  final_referenceindex
    File?  final_referencedict
    File?  final_referencetgz

    String samtools_image         = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String bwa_image              = "biocontainers/bwa:v0.7.17_cv1"
    String gatk_image             = "broadinstitute/gatk3:3.7-0"
    String cov_image              = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/coverage:20230630"
  }

  Boolean realign_bam = defined(final_referencefasta)

  if (realign_bam) {
    call ExtractRGHeaders {
      input:
          inputcram       = inputcram
        , docker          = samtools_image
    }

    scatter (rgheader in ExtractRGHeaders.rgheaders) {
      String rgid = sub(sub(rgheader, "^@RG\\tID:", ""), "\\t.*$", "")
      call ExtractReads {
        input:
            inputcram       = inputcram
          , rgid            = rgid
          , referencefasta  = referencefasta
          , referenceindex  = referenceindex
          , referencedict   = referencedict
          , docker          = samtools_image
      }

      call AlignReads as AlignUnpairedReads {
        input:
            reads          = ExtractReads.reads0
          , rgid           = rgid
          , rgheader       = rgheader
          , referencefasta = select_first([final_referencefasta])
          , referencetgz   = select_first([final_referencetgz])
          , docker         = bwa_image
      }

      call AlignReads as AlignPairedReads {
        input:
            reads          = ExtractReads.reads1
          , reads2         = ExtractReads.reads2
          , rgid           = rgid
          , rgheader       = rgheader
          , referencefasta = select_first([final_referencefasta])
          , referencetgz   = select_first([final_referencetgz])
          , docker         = bwa_image
      }
    }

    scatter (sam in flatten([AlignUnpairedReads.sam, AlignPairedReads.sam])) {
      call SamToBam {
        input:
            sam    = sam
          , docker = samtools_image
      }
    }

    call AddNumbers as GetTotalSize {
      input:
          numbers = SamToBam.size
    }

    call MergeBams {
      input:
          bams    = SamToBam.bam
        , storage = 20 + ceil(3 * GetTotalSize.total)
        , docker = samtools_image
    }
  }

  if (!realign_bam) {
    call CramToBam {
      input:
        input_cram = inputcram,
        sample_name = basename(select_first([inputcram]), ".cram"),
        ref_fasta = referencefasta,
        ref_fasta_index = referenceindex,
        ref_dict = referencedict,
        docker = samtools_image,
        samtools_path = "samtools"
    }
  }

  # call DepthOfCoverageROITask {
  #   input:
  #     ref_fasta = select_first([final_referencefasta, referencefasta]),
  #     ref_fasta_index = select_first([final_referenceindex, referenceindex]),
  #     ref_dict = select_first([final_referencedict, referencedict]),
  #     # bam = select_first([CramToBam.output_bam, SamToBam.realigned_bam]),
  #     # bai = select_first([CramToBam.output_bai, SamToBam.realigned_bai]),
  #     bam = select_first([CramToBam.output_bam, MergeBams.bam]),
  #     bai = select_first([CramToBam.output_bai, MergeBams.bai]),
  #     bed = coveragebed,
  #     max_heap_gb = 31,
  #     docker_image = gatk_image
  # }

  call DepthOfCoverage.DepthOfCoverageWorkflow {
    input:
        run_wgs           = false
      , output_basename   = basename(inputcram, ".cram") + ".cov"
      , ref_fasta         = select_first([final_referencefasta, referencefasta])
      , ref_fasta_index   = select_first([final_referenceindex, referenceindex])
      , ref_dict          = select_first([final_referencedict, referencedict])
      , bam               = select_first([CramToBam.output_bam, MergeBams.bam])
      , bai               = select_first([CramToBam.output_bai, MergeBams.bai])
      , roi_all_bed       = coveragebed
      , roi_genes         = roigenes
      , gene_names        = "gs://lmm-reference-data/roi/HGNC_genenames_05272022.txt"
      , cov_docker_image  = cov_image
      , gatk_docker_image = gatk_image
  }

  output {
    File  roi_sample_interval_summary                = select_first([DepthOfCoverageWorkflow.roi_sample_interval_summary])
    File  roi_sample_summary                         = select_first([DepthOfCoverageWorkflow.roi_sample_summary])
    File? mt_summary                                 = DepthOfCoverageWorkflow.mt_summary
    File  gene_summary                               = select_first([DepthOfCoverageWorkflow.gene_summary])
    File  gene_summary_entrez                        = select_first([DepthOfCoverageWorkflow.gene_summary_entrez])
    File  gene_summary_unknown                       = select_first([DepthOfCoverageWorkflow.gene_summary_unknown])
    File  roi_sample_interval_statistics             = select_first([DepthOfCoverageWorkflow.roi_sample_interval_statistics])
    File  roi_sample_statistics                      = select_first([DepthOfCoverageWorkflow.roi_sample_statistics])
    File  roi_sample_cumulative_coverage_counts      = select_first([DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_counts])
    File  roi_sample_cumulative_coverage_proportions = select_first([DepthOfCoverageWorkflow.roi_sample_cumulative_coverage_proportions])
  }

}

task ExtractRGHeaders {
  input {
    File   inputcram
    String docker
    Int    preemptible = 2
  }

  String OUTPUT  = "rgheaders.txt"
  Int    storage = 20 + ceil(size(inputcram, "GB"))

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  samtools view -H '~{inputcram}' | grep '^@RG' > '~{OUTPUT}'

  INPUTFILES=( '~{inputcram}' )
  ( printf -- '%s\n' "${INPUTFILES[@]}"; find . -type f ) | xargs stat --format=$'%s\t%n'
  ls -lh "${INPUTFILES[@]}"
  ls -lh
  >>>

  output {
    Array[String] rgheaders = read_lines(OUTPUT)
  }

  runtime {
    preemptible : preemptible
    docker      : docker
    disks       : "local-disk ~{storage} HDD"
  }
}

task ExtractReads {
  input {
    File   inputcram
    String rgid
    File   referencefasta
    File   referenceindex
    File   referencedict
    Int    compressionlevel = 6
    String docker
    Int    preemptible      = 2
    Int    memory           = 16
  }

  Float referencesize = size(referencefasta, "GB") +
                        size(referenceindex, "GB") +
                        size(referencedict,  "GB")
  # Int   storage       = 20 + 3 * ceil(size(inputcram, "GB") + referencesize)
  Int   storage       = 100

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  REFERENCEFASTA='~{referencefasta}'

  if [[ ~{referenceindex} != ${REFERENCEFASTA}.fai ]]
  then
    ln --symbolic '~{referenceindex}' "${REFERENCEFASTA}.fai"
  fi

  if [[ ~{referencedict} != ${REFERENCEFASTA%.*}.dict ]]
  then
    ln --symbolic '~{referencedict}' "${REFERENCEFASTA%.*}.dict"
  fi

    samtools                            \
        view                            \
        -h                              \
        --reference "${REFERENCEFASTA}" \
        '~{inputcram}'                  \
  | perl -lne '
      if ( /^\@/ ) { print; next; }
      print if /\tRG:Z:\Q~{rgid}\E$/;
    '                                   \
  | samtools                            \
        fastq                           \
        -c ~{compressionlevel}          \
        --reference "${REFERENCEFASTA}" \
        -1 reads1.fastq.gz              \
        -2 reads2.fastq.gz              \
        -s reads0.fastq.gz              \
        -

  INPUTFILES=( '~{inputcram}' "${REFERENCEFASTA}" '~{referenceindex}' '~{referencedict}' )
  ( printf -- '%s\n' "${INPUTFILES[@]}"; find . -type f ) | xargs stat --dereference --format=$'%s\t%n'
  ls -lh "${INPUTFILES[@]}"
  ls -lh
  >>>

  output {
    File reads0 = "reads0.fastq.gz"
    File reads1 = "reads1.fastq.gz"
    File reads2 = "reads2.fastq.gz"
  }

  runtime {
    preemptible : preemptible
    docker      : docker
    memory      : "~{memory} GB"
    disks       : "local-disk ~{storage} HDD"
  }
}

task AlignReads {
  input {
    File   reads
    File?  reads2
    String rgid
    String rgheader
    File   referencefasta
    File   referencetgz
    String docker
    Int    memory         = 16
    Int    ncpus          = 96
    Int    preemptible    = 2
  }

  Float  readssize     = size(reads, "GB") +
                         if defined(reads2) then size(reads2, "GB") else 0

  Float  referencesize = size(referencefasta, "GB") +
                         size(referencetgz  , "GB") * 2.5

  Int    storage       = 20 + ceil(2 * readssize + referencesize)

  String suffix        = if defined(reads2) then "paired" else "unpaired"
  String SAM           = "~{rgid}__~{suffix}.sam"
  command <<<

  memory_used() {
    free | head --lines=2 | tail --lines=1 | tr --squeeze ' ' $'\t' | cut --fields=3
  }

  main() {
    (
      # -----------------------------------------------------------------

      bwa index >&2
      bwa mem   >&2

      # -----------------------------------------------------------------

      set -o errexit
      set -o pipefail
      set -o nounset
      set -o xtrace

      # -----------------------------------------------------------------
      # The following is required because `bwa index` needs to be able to
      # write files to its argument parent directory.
      mkdir reference
      REFERENCEFASTA="reference/$( basename '~{referencefasta}' )"
      ln --symbolic '~{referencefasta}' "${REFERENCEFASTA}"
      # -----------------------------------------------------------------

      # bwa index -a bwtsw "${REFERENCEFASTA}"

      (
        cd reference
        tar xvzf '~{referencetgz}'
      )

      if ~{if defined(reads2) then "true" else "false"}
      then
          READS=( '~{reads}' '~{reads2}' )
      else
          READS=( '~{reads}' )
      fi

      RGHEADER="$( perl -lpe 's/\t/\\t/g' <<<'~{rgheader}' )"

      NPROC="$( nproc )"
      bwa                     \
          mem                 \
          -K 100000000        \
          -t "${NPROC}"       \
          -R "${RGHEADER}"    \
          "${REFERENCEFASTA}" \
          "${READS[@]}"       \
        > '~{SAM}'

      INPUTFILES=( "${REFERENCEFASTA}" "${READS[@]}" )
      ( printf -- '%s\n' "${INPUTFILES[@]}"; find . \! -type d ) | xargs stat --dereference --format=$'%s\t%n'
      ls -lh "${INPUTFILES[@]}"
      ls -lh
    )
  }

  ( main ) &
  mainpid=$!

  MINIMUM_SLEEPTIME=1
  MAXIMUM_SLEEPTIME=256
  SLEEPTIME="${MINIMUM_SLEEPTIME}"
  MAXIMUM_MEMORY_USED=0
  while true
  do
      MEMORY_USED="$( memory_used )"
      if (( MAXIMUM_MEMORY_USED < MEMORY_USED))
      then
          free --human >&2
          MAXIMUM_MEMORY_USED="${MEMORY_USED}"
          if (( SLEEPTIME > MINIMUM_SLEEPTIME ))
          then
              SLEEPTIME=$(( SLEEPTIME / 2 ))
          fi
      elif (( SLEEPTIME < MAXIMUM_SLEEPTIME ))
      then
          SLEEPTIME=$(( SLEEPTIME * 2 ))
      fi
      if ! kill -0 "${mainpid}" 2>/dev/null
      then
          wait "${mainpid}"
          exit_status=$?
          printf -- 'main exited with status %d\n' "${exit_status}" >&2
          printf -- 'maximum memory used: %d\n' "${MAXIMUM_MEMORY_USED}" >&2
          exit "${exit_status}"
      fi
      sleep "${SLEEPTIME}"
  done
  >>>

  output {
    File sam = SAM
  }

  runtime {
    preemptible : preemptible
    docker      : docker
    memory      : "~{memory} GB"
    disks       : "local-disk ~{storage} HDD"
    cpu         : ncpus
  }
}

task SamToBam {
  input {
    File   sam
    String docker
    Int    preemptible = 2
    Int    memory      = 16
  }

  Float  samsize   = size(sam, "GB")
  Int    storage   = 20 + ceil(4 * samsize)
  String BAM       = basename(sam, ".sam") + ".bam"

  command <<<
  df
  df --human

  printf -- '\n\n'
  printf -- 'samsize: %5.1f GB\n' '~{samsize}'
  printf -- 'storage: %3d   GB\n' '~{storage}'
  printf -- '\n\n'

  # ---------------------------------------------------------------

  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  samtools view -Sb '~{sam}' > '~{BAM}'

  # -----------------------------------------------------------------
  df
  df --human
  printf -- '\n\n'

  INPUTFILES=( '~{sam}' )
  ( printf -- '%s\n' "${INPUTFILES[@]}"; find . \! -type d ) | xargs stat --format=$'%s\t%n'
  ls -lh "${INPUTFILES[@]}"
  ls -lh
  >>>

  output {
    File bam  = BAM
    Int  size = ceil(size(BAM, "GB"))
  }

  runtime {
    preemptible : preemptible
    docker      : docker
    memory      : "~{memory} GB"
    disks       : "local-disk ~{storage} HDD"
  }
}


task AddNumbers {
  input {
    Array[Float] numbers
  }

  String OUTPUT = "OUTPUT.txt"
  command <<<
  python3 <<EOF
  numbers = [~{sep=", " numbers}]
  with open('~{OUTPUT}', 'w') as writer:
      print(str(sum(numbers)), file=writer)
  EOF
  >>>

  output {
    Float total = read_float(OUTPUT)
  }

  runtime {
    docker: "python:3.11"
  }
}


task MergeBams {
  input {
    Array[File] bams
    Int         storage
    String      docker
    Int         preemptible = 2
    Int         memory      = 16
  }

  String BAM = "realigned.bam"
  String BAI = "realigned.bai"

  command <<<
  df
  df --human

  printf -- '\n\n'
  printf -- 'storage: %3d   GB\n' '~{storage}'
  printf -- '\n\n'

  # ---------------------------------------------------------------

  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  samtools merge unsorted.bam '~{sep="' '" bams}'
  samtools sort -o '~{BAM}' unsorted.bam
  samtools index '~{BAM}' '~{BAI}'

  # -----------------------------------------------------------------
  df
  df --human
  printf -- '\n\n'

  INPUTFILES=( '~{sep="' '" bams}' )
  ( printf -- '%s\n' "${INPUTFILES[@]}"; find . \! -type d ) | xargs stat --format=$'%s\t%n'
  ls -lh "${INPUTFILES[@]}"
  ls -lh
  >>>

  output {
    File bam = BAM
    File bai = BAI
  }

  runtime {
    preemptible : preemptible
    docker      : docker
    memory      : "~{memory} GB"
    disks       : "local-disk ~{storage} HDD"
  }
}



task CramToBam {
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

# task DepthOfCoverageROITask {
#     input {
#         File ref_fasta
#         File ref_fasta_index
#         File ref_dict
#         File bam
#         File bai
#         File bed
#         String output_basename = sub(basename(bam), "\\.(bam|BAM|cram|CRAM)$", "") + ".cov"
#         Int max_heap_gb
#         Int disk_size = ceil(size(bam, "GB") * 1.5) + 10
#         String docker_image
#         Int preemptible = 1
#     }
#
#     command <<<
#         set -euxo pipefail
#         mkdir cov_out
#         java -Xmx~{max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
#             -I "~{bam}" \
#             -ct 8 -ct 15 \
#             -R "~{ref_fasta}" \
#             -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
#             --printBaseCounts \
#             -o "cov_out/~{output_basename}.roibed" \
#             -L "~{bed}"
#         ls -l cov_out
#     >>>
#
#     runtime {
#         docker: "~{docker_image}"
#         memory: (max_heap_gb + 4) + "GB"
#         disks: "local-disk " + disk_size + " SSD"
#         preemptible: preemptible
#     }
#
#     output {
#         File sample_interval_summary = "cov_out/~{output_basename}.roibed.sample_interval_summary"
#         File sample_interval_statistics = "cov_out/~{output_basename}.roibed.sample_interval_statistics"
#         File sample_statistics = "cov_out/~{output_basename}.roibed.sample_statistics"
#         File sample_summary = "cov_out/~{output_basename}.roibed.sample_summary"
#         File sample_cumulative_coverage_counts = "cov_out/~{output_basename}.roibed.sample_cumulative_coverage_counts"
#         File sample_cumulative_coverage_proportions = "cov_out/~{output_basename}.roibed.sample_cumulative_coverage_proportions"
#     }
# }

  # call CramToBam {
  #   input:
  #       inputcram      = inputcram
  #     , referencedict  = referencedict
  #     , referencefasta = referencefasta
  #     , referenceindex   = referenceindex
  #     , docker         = samtools_image
  #     , samtools_path  = "samtools"
  # }

  # call DepthOfCoverageROITask {
  #   input:
  #       bam            = CramToBam.bam
  #     , bai            = CramToBam.bai
  #     , bed            = coveragebed
  #     , referencedict  = referencedict
  #     , referencefasta = referencefasta
  #     , referenceindex   = referenceindex

  #     , max_heap_gb    = 31
  #     , docker         = gatk_image
  # }

# task CramToBam {
#   input {
#     # Command parameters
#     File   inputcram
#     File   referencefasta
#     File   referenceindex
#     File   referencedict

#     # Runtime parameters
#     String  docker
#     Int?    memory
#     Int?    storage
#     Int?    preemptible
#     String  samtools_path
#   }

#   String sample_name    = basename(inputcram, ".cram")
#   Float output_bam_size = size(inputcram, "GB") / 0.40
#   Float referencesize   = size(referencefasta, "GB") +
#                           size(referenceindex  , "GB") +
#                           size(referencedict , "GB")
#   Int   disk_size       = 20 + ceil(size(inputcram, "GB") +
#                                     output_bam_size +
#                                     referencesize)

#   command {
#     set -e
#     set -o pipefail

#     ~{samtools_path}                      \
#         view                              \
#         --with-header                     \
#         --reference '~{referencefasta}'   \
#         '~{inputcram}'                    \
#       | ~{samtools_path}                  \
#             view                          \
#             --bam                         \
#             --output '~{sample_name}.bam' \
#             -

#     ~{samtools_path}                  \
#         index                         \
#         --bai                         \
#         --output '~{sample_name}.bai' \
#         '~{sample_name}.bam'
#   }
#   runtime {
#     docker: docker
#     memory: "~{select_first([memory, 16])} GB"
#     disks: "local-disk ~{select_first([storage, disk_size])} HDD"
#     preemptible: select_first([preemptible, 2])
#   }

#   output {
#     File bam = "~{sample_name}.bam"
#     File bai = "~{sample_name}.bai"
#   }
# }

# task DepthOfCoverageROITask {
#     input {
#         File   bam
#         File   bai
#         File   bed
#         File   referencefasta
#         File   referenceindex
#         File   referencedict
#         String output_basename = sub(basename(bam), "\\.(bam|BAM|cram|CRAM)$", "") + ".cov"
#         Int    max_heap_gb
#         Int    disk_size = ceil(size(bam, "GB") * 1.5) + 10
#         String docker
#         Int    preemptible = 1
#     }

#     command <<<
#     set -euxo pipefail

#     mkdir cov_out

#     java                                              \
#         -Xmx~{max_heap_gb}g                           \
#         -jar /usr/GenomeAnalysisTK.jar                \
#         -T DepthOfCoverage                            \
#         -I "~{bam}"                                   \
#         -ct 8 -ct 15                                  \
#         -R "~{referencefasta}"                        \
#         -dt BY_SAMPLE                                 \
#         -dcov 1000                                    \
#         -l INFO                                       \
#         --omitDepthOutputAtEachBase                   \
#         --minBaseQuality 10                           \
#         --minMappingQuality 17                        \
#         --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
#         --printBaseCounts                             \
#         -o "cov_out/~{output_basename}.roibed"        \
#         -L "~{bed}"

#     ls -l cov_out
#     >>>

#     runtime {
#         docker      : docker
#         memory      : "~{max_heap_gb + 4} GB"
#         disks       : "local-disk {disk_size} SSD"
#         preemptible : preemptible
#     }

#     output {
#         File sample_interval_summary = "cov_out/~{output_basename}.roibed.sample_interval_summary"
#         File sample_interval_statistics = "cov_out/~{output_basename}.roibed.sample_interval_statistics"
#         File sample_statistics = "cov_out/~{output_basename}.roibed.sample_statistics"
#         File sample_summary = "cov_out/~{output_basename}.roibed.sample_summary"
#         File sample_cumulative_coverage_counts = "cov_out/~{output_basename}.roibed.sample_cumulative_coverage_counts"
#         File sample_cumulative_coverage_proportions = "cov_out/~{output_basename}.roibed.sample_cumulative_coverage_proportions"
#     }
# }
