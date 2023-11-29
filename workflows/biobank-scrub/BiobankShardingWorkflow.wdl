version 1.0

# -----------------------------------------------------------------------------

workflow BiobankShardingWorkflow {
  input {
  # Int      shard_area       = 50000000
    Int      shard_area       = 500000000
    # File     areas_tsv
    String   areas_tsv
    String   base_directory
    String   docker_image     = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/biobank-scrub:1.0.1"

    # -------------------------------------------------------------------------
    # NB: The parameters below this point are rarely needed!

    Int      do_excerpts_code = 0  # SHOULD NEVER BE USED IN PRODUCTION; used
                                   # only when generating fixtures for testing

    String   output_base      = "" # If do_excerpts_code > 0, output_base
                                   # MUST be specified as a non-empty string
                                   # different from base_directory.
  }

  # call ShowEnvironment {
  #   input:
  #     docker_image = docker_image
  # }

  # ---------------------------------------------------------------------------

  call MakeShardingBatches {
    input:
      shard_area     = shard_area,
      areas_tsv      = areas_tsv,
      base_directory = base_directory,
      docker_image   = docker_image
  }

  scatter (storage_and_batch in MakeShardingBatches.storage_and_batch_array) {
    call ShardVcfs {
      input:
        base_directory   = base_directory,
        batch            = storage_and_batch.batch,
        do_excerpts_code = do_excerpts_code,
        output_base      = output_base,
        storage          = storage_and_batch.storage,
        docker_image     = docker_image
    }
  }

  call ConsolidateBatchedResults {
    input:
      batched_results_array = ShardVcfs.results,
      docker_image          = docker_image
  }

  output {
    # File   environment      = ShowEnvironment.environment
    File   sharding_results = ConsolidateBatchedResults.results
  }
}

# -----------------------------------------------------------------------------

# task ShowEnvironment {
#   input {
#     String  docker_image
#   }
# 
#   String   OUTPUTDIR = "OUTPUT"
#   String   STDOUT    = OUTPUTDIR + "/STDOUT"
# 
#   command <<<
#   set -o errexit
#   set -o pipefail
#   # set -o nounset
#   set -o xtrace
#   # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
# 
#   show_environment() {
#       typeset -p
#       rclone version
#   }
# 
#   export RCLONE_LOG_LEVEL=DEBUG
#   export RCLONE_STATS_LOG_LEVEL=DEBUG
# 
#   ### FIXME: DELETE THE FOLLOWING LINE!!!
#   /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w vcf-merge-workflow-test -r gs://fc-secure-209c8de0-8bee-473e-bda5-2b0e8f1691b1
# 
#   mkdir --parents '~{OUTPUTDIR}'
# 
#   exec >&2
#   show_environment | tee '~{STDOUT}'
# 
#   >>>
# 
#   output {
#     File environment = STDOUT
#   }
# 
#   runtime {
#     docker: docker_image
#   }
# }


task  MakeShardingBatches {
  input {
    Int      shard_area
    # File     areas_tsv
    String   areas_tsv
    String   base_directory
    String   docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{areas_tsv}'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{base_directory}'

  AREAS_TSV=/tmp/areas.tsv

  # rclone copyto '~{areas_tsv}' "${AREAS_TSV}"
  gsutil cp '~{areas_tsv}' "${AREAS_TSV}"

  mkdir --parents '~{OUTPUTDIR}'

  # make_sharding_batches.py \
  #     '~{shard_area}'      \
  #     '~{areas_tsv}'       \
  #     '~{base_directory}'  \
  #   > '~{STDOUT}'

  make_sharding_batches.py \
      '~{shard_area}'      \
      "${AREAS_TSV}"       \
      '~{base_directory}'  \
    > '~{STDOUT}'


  >>>

  output {
    Array[Object]   storage_and_batch_array = read_json(STDOUT)
  }

  runtime {
    docker: docker_image
  }
}


task  ShardVcfs {
  # ...
  input {
    String          base_directory
    Array[Object]   batch
    Int             do_excerpts_code
    String          output_base
    Int             storage
    String          docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  export RCLONE_LOG_LEVEL=DEBUG
  export RCLONE_STATS_LOG_LEVEL=DEBUG

  (
    exec >&2
    typeset -p
    pwd
    ps -e --format=pid,ppid,%cpu,%mem,vsz,rss,tt,stat,start,etime,time,cmd
    pstree --arguments --compact-not --long --thread-names --show-parents --show-pids --numeric-sort --ascii
  )

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{base_directory}'

  mkdir --parents '~{OUTPUTDIR}'

  declare -a EXTRAFLAGS
  EXTRAFLAGS=()

  DO_EXCERPTS_CODE=~{do_excerpts_code}
  OUTPUT_BASE='~{output_base}'

  if (( DO_EXCERPTS_CODE != 0 ))
  then
      EXTRAFLAGS+=( "--do-excerpts-code=${DO_EXCERPTS_CODE}" )
  fi

  if [[ -n ${OUTPUT_BASE} ]]
  then
      EXTRAFLAGS+=( "--output-base=${OUTPUT_BASE}" )
  elif (( DO_EXCERPTS_CODE > 0 ))
  then
      printf --                                                         \
          'ERROR: positive DO_EXCERPTS_CODE (%d) and empty OUTPUT_BASE' \
          "${DO_EXCERPTS_CODE}" >&2
      exit 1
  fi

  shard_vcfs.py                           \
      --logging-level=DEBUG               \
      ${EXTRAFLAGS[@]+"${EXTRAFLAGS[@]}"} \
      '~{base_directory}'                 \
      '~{write_json(batch)}'              \
  > '~{STDOUT}'

  >>>

  output {
    # String   sequencing_dummy = read_string(DUMMY)
    Object   results          = read_json(STDOUT)
  }

  runtime {
    docker: docker_image
    disks:  "local-disk ~{storage} HDD"
  }
}


task  ConsolidateBatchedResults {
  # ...
  input {
    Array[Object]   batched_results_array
    String          docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  mkdir --parents '~{OUTPUTDIR}'

  consolidate_batched_results.py             \
      '~{write_json(batched_results_array)}' \
    > '~{STDOUT}'

  >>>

  output {
    File   results = STDOUT
  }

  runtime {
    docker: docker_image
  }
}

# -----------------------------------------------------------------------------