version 1.0

import "https://raw.githubusercontent.com/mgbpm/biofx-workflows/refs/heads/main/steps/Utilities.wdl"

workflow T1kHlaTyping {
  input {
    File    bam_or_cram
    File    bai_or_crai
    File?   reference_fasta
    File?   reference_fai
    String  preset
    Int?    memorygb
    Int?    storagegb
    String  docker          = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/t1k:0.0.1"
  }

  File default_reference_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File default_reference_fai   = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  String format = if      (sub(bam_or_cram, "^.*\\.[Bb][Aa][Mm]$",     "") == "") then "BAM" 
                  else if (sub(bam_or_cram, "^.*\\.[Cc][Rr][Aa][Mm]$", "") == "") then "CRAM"
                  else "UNSUPPORTED"

  if (format == "UNSUPPORTED") {
    call Utilities.FailTask as UnsupportedFormatError {
      input:
        error_message = "the extension of ~{bam_or_cram} is neither .bam nor .cram"
    }
  }

  if (format == "BAM" || format == "CRAM") {

    if (format == "CRAM") {
      Boolean reference_provided = defined(reference_fasta) && defined(reference_fai)
      File maybe_reference_fasta = if reference_provided then select_first([reference_fasta]) else default_reference_fasta
      File maybe_reference_fai   = if reference_provided then select_first([reference_fai])   else default_reference_fai
    }

    # Float storage_multiplier = if format == "CRAM" then 3.5 else 1.5
    Int storage_multiplier = if format == "CRAM" then 4 else 2

    call T1kHlaTypingTask {
      input:
          format          = format
        , bam_or_cram     = bam_or_cram
        , bai_or_crai     = bai_or_crai
        , reference_fasta = maybe_reference_fasta
        , reference_fai   = maybe_reference_fai
        , preset          = preset
        , memorygb        = select_first([memorygb,  ceil(16 + size(bam_or_cram, "GB") * 0.5)])
        , storagegb       = select_first([storagegb, ceil( 8 + size(bam_or_cram, "GB") * storage_multiplier
                                                             + size(maybe_reference_fasta, "GB"))])
        , docker          = docker
    }
  }

  output {
    Array[File] outputfiles = select_first([T1kHlaTypingTask.outputfiles, []])
  }
}

task T1kHlaTypingTask {
  input {
    String  format
    File    bam_or_cram
    File    bai_or_crai
    File?   reference_fasta
    File?   reference_fai
    String  preset
    String  basename_pattern     = "\\.[^\\.]+$"
    String  basename_replacement = ""
    String  docker
    Int     memorygb
    Int     storagegb
  }

  String OUTPUT = "OUTPUT"
  String STEM   = sub(basename(bam_or_cram),
                      basename_pattern,
                      basename_replacement)
  String PREFIX = "~{OUTPUT}/~{STEM}"

  command <<<
  footprint () {
      (
          set +o errexit

          exec >&2 2>/dev/null

          du --summarize                 --block-size=1  "${PWD}"
          du --summarize --apparent-size --block-size=G  "${PWD}"
          du --summarize                 --block-size=G  "${PWD}"
          du --summarize --apparent-size --block-size=GB "${PWD}"
          du --summarize                 --block-size=GB "${PWD}"
          du --summarize --apparent-size --human         "${PWD}"
          du --summarize                 --human         "${PWD}"

          printf -- '\n'
          df         "${PWD}"
          df --human "${PWD}"
      )

      true
  }

  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  (
      exec >&2
      printf -- 'memorygb : %s\n' '~{memorygb}'
      printf -- 'storagegb: %s\n' '~{storagegb}'
      printf -- '\n'
  )

  footprint

  # --------------------------------------------------------------------------

  rm --recursive --force '~{OUTPUT}'
  mkdir --parents '~{OUTPUT}'

  TEMPORARY="$( mktemp --directory )"

  INPUTBAM="${TEMPORARY}/input.bam"
  INPUTBAI="${TEMPORARY}/input.bai"

  if [[ '~{format}' == 'CRAM' ]]
  then
      WANTEDMAGIC=CRAM
      WANTEDTYPE=cram

      INPUTCRAM="${TEMPORARY}/input.cram"
      ln --symbolic '~{bam_or_cram}' "${INPUTCRAM}"
      ln --symbolic '~{bai_or_crai}' "${TEMPORARY}/input.crai"
  else
      WANTEDMAGIC=$'BAM\001'
      WANTEDTYPE=bam

      ln --symbolic '~{bam_or_cram}' "${INPUTBAM}"
      ln --symbolic '~{bai_or_crai}' "${INPUTBAI}"
  fi

  if [[ "${WANTEDMAGIC}" != "$( zcat --force '~{bam_or_cram}' | head --bytes 4 )" ]]
  then
      printf -- 'ERROR: %s is not a %s file\n' '~{bam_or_cram}' "${WANTEDTYPE}" >&2
      exit 1
  fi

  if [[ '~{format}' == 'CRAM' ]]
  then
      REFERENCE_FASTA="${TEMPORARY}/reference.fasta"
      REFERENCE_FAI="${TEMPORARY}/reference.fai"
      ln --symbolic '~{reference_fasta}' "${REFERENCE_FASTA}"
      ln --symbolic '~{reference_fai}'   "${REFERENCE_FAI}"

      NTHREADS="$( nproc )"
      /usr/bin/time -v                        \
      samtools                                \
          view                                \
          --bam                               \
          --threads      "${NTHREADS}"        \
          --reference    "${REFERENCE_FASTA}" \
          --output       "${INPUTBAM}"        \
          "${INPUTCRAM}"

      /usr/bin/time -v                \
      samtools                        \
          index                       \
          --bai                       \
          --threads     "${NTHREADS}" \
          --output      "${INPUTBAI}" \
          "${INPUTBAM}"
  fi

  # --------------------------------------------------------------------------

  mkdir --parents "$( dirname '~{PREFIX}' )"

  /usr/bin/time -v t1k '~{preset}' "${INPUTBAM}" '~{PREFIX}'

  # ----------------------------------------------------------------------------

  footprint

  >>>

  output {
    Array[File] outputfiles = glob("~{PREFIX}*")
  }

  runtime {
    disks:  "local-disk ~{storagegb} SSD"
    memory: "~{memorygb}GB"
    docker: docker
  }
}
