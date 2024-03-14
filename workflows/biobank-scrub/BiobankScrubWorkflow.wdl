version 1.0
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# -----------------------------------------------------------------------------

import "BiobankUtils.wdl"

# -----------------------------------------------------------------------------

struct RundirInitialized {
  String          staging_area
  File?           maybe_saved_withdrawn_list
  Int             attempt_index
}

# -----------------------------------------------------------------------------

workflow BiobankScrubWorkflow {
  input {
    String   runid
    Int      nbatches                = 500
    String   scrubsdir               = "gs://mgbpm-biobank-data/scrubs"
    String   initial_datadir         = "gs://mgbpm-biobank-data/datasets/initial"
    String   current_datadir         = "gs://mgbpm-biobank-data/datasets/current"
    String?  release_datadir
    File     latest_withdrawn_list   = "gs://mgbpm-biobank-data/reference/biobank/withdrawn.txt"
    File?    withdrawn_list_override
    Boolean  force                   = false
    Int      scrub_memory            = 26   # in GB
    Int      scrub_disk_size         = 375  # in GB
    String   docker_image            = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/biobank-scrub:1.0.8"
  }

  call BiobankUtils.ShowEnvironment as ShowEnvironment {
    input:
      docker_image = docker_image
  }

  Object inputs = object {
    runid                   : runid,
    nbatches                : nbatches,
    scrubsdir               : scrubsdir,
    initial_datadir         : initial_datadir,
    current_datadir         : current_datadir,
    release_datadir         : release_datadir,
    latest_withdrawn_list   : latest_withdrawn_list,
    withdrawn_list_override : withdrawn_list_override,
    force                   : force,
    scrub_memory            : scrub_memory,
    scrub_disk_size         : scrub_disk_size,
    docker_image            : docker_image
  }

  call ValidateInputs {
    input:
      inputs = inputs
  }

  String rundir = scrubsdir + '/' + runid

  call MaybeInitializeRundir {
    input:
      rundir           = rundir,
      inputs           = inputs,
      sequencing_dummy = ValidateInputs.sequencing_dummy
  }

  String   staging_area               = MaybeInitializeRundir.info.staging_area
  Int      attempt_index              = MaybeInitializeRundir.info.attempt_index
  File?    maybe_saved_withdrawn_list =
    MaybeInitializeRundir.info.maybe_saved_withdrawn_list

  File withdrawn_list = select_first([
    withdrawn_list_override,
    maybe_saved_withdrawn_list,
    latest_withdrawn_list
  ])

  call ListDatasetIds {
    input:
      base_datadir = current_datadir,
      docker_image = docker_image
  }

  scatter (dataset_id in ListDatasetIds.dataset_ids) {
    call FindNonCompliant {
      input:
        dataset_id      = dataset_id,
        withdrawn_list  = withdrawn_list,
        current_datadir = current_datadir,
        initial_datadir = initial_datadir,
        docker_image    = docker_image
    }
  }

  Array[Object] non_compliant  = flatten(FindNonCompliant.non_compliant)

  call MakeScrubBatches {
    input:
      non_compliant = non_compliant,
      nbatches      = nbatches,
      rundir        = rundir,         # to check for already-scrubbed files
      docker_image  = docker_image
  }

  scatter (batch in MakeScrubBatches.batches) {
    call ScrubBatch {
      input:
        withdrawn_list  = withdrawn_list,
        initial_datadir = initial_datadir,
        rundir          = rundir,
        batch           = batch,
        memory          = "~{scrub_memory}GB",
        storage         = "local-disk ~{scrub_disk_size} HDD",
        docker_image    = docker_image
    }
  }

  # -------------------------------------------------------------------

  call CollectShards {
    input:
      rundir              = rundir,
      non_compliant       = non_compliant,
      sequencing_sentinel = ScrubBatch.results,  # ignored
      docker_image        = docker_image
  }

  scatter (storage_and_batch in CollectShards.storage_and_batch_array) {
    call ConcatenateShards {
      input:
        rundir       = rundir,
        batch        = storage_and_batch.batch,
        storage      = storage_and_batch.storage,
        docker_image = docker_image
    }
  }

  call CheckBatchedResults as CheckScrubResults {
    input:
      batched_results = ScrubBatch.results,
      docker_image    = docker_image
  }

  call CheckBatchedResults as CheckConcatenationResults {
    input:
      batched_results = ConcatenateShards.results,
      docker_image    = docker_image
  }

  if (CheckScrubResults.isok && CheckConcatenationResults.isok) {

    call MakePushBatches {
      input:
        rundir       = rundir,
        docker_image = docker_image
    }

    scatter (batch in MakePushBatches.batches) {
      call PushScrubbed {
        input:
          rundir          = rundir,
          release_datadir = select_first([
                                          release_datadir,
                                          current_datadir
                                         ]),
          batch           = batch,
          docker_image    = docker_image
      }
    }

  }

  call Summarize {
    input:
      scrub_results         = ScrubBatch.results,
      concatenation_results = ConcatenateShards.results,
      push_results          = select_first([PushScrubbed.results, []]),
      docker_image          = docker_image
  }

  output {
    File   environment = ShowEnvironment.environment
    File   summary     = Summarize.results
  }
}

# -----------------------------------------------------------------------------

task  ValidateInputs {
  input {
    Object inputs
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"
  String   DUMMY     = OUTPUTDIR + "/DUMMY"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  mkdir --parents '~{OUTPUTDIR}'

  ( hostname --long ; date ) > '~{DUMMY}'

  validate_inputs.py          \
      --logging-level=DEBUG   \
      '~{write_json(inputs)}' \
    | tee '~{STDOUT}'

  date >> '~{DUMMY}'
  >>>

  output {
    String   sequencing_dummy = read_string(DUMMY)
  }

  runtime {
    preemptible: 2
    docker:      inputs.docker_image
  }
}

task  MaybeInitializeRundir {
  input {
    String rundir
    Object inputs
    String sequencing_dummy
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  export RCLONE_LOG_LEVEL=DEBUG
  export RCLONE_STATS_LOG_LEVEL=DEBUG

  typeset -p >&2
  rclone version >&2

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  mkdir --parents '~{OUTPUTDIR}'

  maybe_initialize_rundir.py  \
      --logging-level=DEBUG   \
      '~{rundir}'             \
      '~{write_json(inputs)}' \
      | tee '~{STDOUT}'
  >>>

  output {
    RundirInitialized info = read_json(STDOUT)
  }

  runtime {
    # NO PREEMPTION FOR THIS TASK
    preemptible: 0
    docker:      inputs.docker_image
  }
}


task  ListDatasetIds {
  input {
    String  base_datadir
    String  docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{base_datadir}'

  mkdir --parents '~{OUTPUTDIR}'
  list_dataset_ids.py '~{base_datadir}' | tee '~{STDOUT}'
  >>>

  output {
    Array[String]   dataset_ids = read_lines(STDOUT)
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}


task  FindNonCompliant {
  input {
    String  current_datadir  # url to current data directory
    String  initial_datadir  # url to initial data directory
    File    withdrawn_list   # file containing all withdrawn subjects,
                             # one-per-line
    String  dataset_id       # 4-digit identifier
    String  docker_image
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

  typeset -p >&2
  rclone version >&2

  # LOCATIONS=( '~{withdrawn_list}' '~{current_datadir}' '~{initial_datadir}' )
  LOCATIONS=( '~{initial_datadir}' )
  export PATH="${PACKAGESDIR}/biofx-orchestration-utils/bin:${PATH}"
  for LOCATION in "${LOCATIONS[@]}"
  do
      /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r "${LOCATION}"
  done

  mkdir --parents '~{OUTPUTDIR}'
  find_non_compliant.py    \
      '~{withdrawn_list}'  \
      '~{current_datadir}' \
      '~{initial_datadir}' \
      '~{dataset_id}'      \
    | tee '~{STDOUT}'
  >>>

  output {
    Array[Object]   non_compliant = read_json(STDOUT)
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}


task  MakeScrubBatches {
  input {
    Array[Object]   non_compliant   # each object in this array has
                                    # members path and type
    Int             nbatches        # desired number of batches
    String          rundir          # url to a run directory
    String          docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  mkdir --parents '~{OUTPUTDIR}'
  make_scrub_batches.py              \
      --logging-level=DEBUG          \
      '~{write_json(non_compliant)}' \
      ~{nbatches}                    \
      '~{rundir}'                    \
    | tee '~{STDOUT}'
  # PRINTS TO STDOUT AN ARRAY OF nbatches ARRAYS OF OBJECTS; EACH OBJECT IN
  # THIS OUTPUT CORRESPONDS TO A (STILL-PENDING) OBJECT IN THE non_compliant
  # INPUT.
  >>>

  output {
    Array[Array[Object]]   batches = read_json(STDOUT)
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}


task  ScrubBatch {
  input {
    File            withdrawn_list
    String          initial_datadir
    String          rundir
    Array[Object]   batch
    String          memory
    String          storage
    String          docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"
  String   DUMMY     = OUTPUTDIR + "/DUMMY"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  export STARTDIR="${PWD}"

  (
    exec >&2

    free --bytes
    free --human
  )

  (
    exec >&2

    df --block-size=1  "${TMPDIR}"
    df --block-size=G  "${TMPDIR}"
    df --block-size=GB "${TMPDIR}"
    df --human         "${TMPDIR}"
  )

  (
    exec >&2

    set +o errexit

    du --summarize --apparent-size --block-size=1  "${STARTDIR}" 2>/dev/null
    du --summarize                 --block-size=1  "${STARTDIR}" 2>/dev/null
    du --summarize --apparent-size --block-size=G  "${STARTDIR}" 2>/dev/null
    du --summarize                 --block-size=G  "${STARTDIR}" 2>/dev/null
    du --summarize --apparent-size --block-size=GB "${STARTDIR}" 2>/dev/null
    du --summarize                 --block-size=GB "${STARTDIR}" 2>/dev/null
    du --summarize --apparent-size --human         "${STARTDIR}" 2>/dev/null
    du --summarize                 --human         "${STARTDIR}" 2>/dev/null

    true
  )

  export RCLONE_LOG_LEVEL=DEBUG
  export RCLONE_STATS_LOG_LEVEL=DEBUG

  (
      exec >&2
      typeset -p
      rclone version
      which scrub_batch.py
      BINDIR="$( dirname "$( which scrub_batch.py )" )"
      ls -Altrd "${BINDIR}"
      ls -Altr  "${BINDIR}"
  )

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  mkdir --parents '~{OUTPUTDIR}'
  ( hostname --long ; date ) > '~{DUMMY}'

  /usr/bin/time -v               \
      scrub_batch.py             \
          --logging-level=DEBUG  \
          '~{withdrawn_list}'    \
          '~{initial_datadir}'   \
          '~{rundir}'            \
          '~{write_json(batch)}' \
        | tee '~{STDOUT}'

  date >> '~{DUMMY}'

  (
    exec >&2

    df --block-size=1  "${TMPDIR}"
    df --block-size=G  "${TMPDIR}"
    df --block-size=GB "${TMPDIR}"
    df --human         "${TMPDIR}"
  )

  (
    exec >&2

    set +o errexit

    du --summarize --apparent-size --block-size=1  "${STARTDIR}" 2>/dev/null
    du --summarize                 --block-size=1  "${STARTDIR}" 2>/dev/null
    du --summarize --apparent-size --block-size=G  "${STARTDIR}" 2>/dev/null
    du --summarize                 --block-size=G  "${STARTDIR}" 2>/dev/null
    du --summarize --apparent-size --block-size=GB "${STARTDIR}" 2>/dev/null
    du --summarize                 --block-size=GB "${STARTDIR}" 2>/dev/null
    du --summarize --apparent-size --human         "${STARTDIR}" 2>/dev/null
    du --summarize                 --human         "${STARTDIR}" 2>/dev/null

    true
  )

  >>>

  output {
    Object   results = read_json(STDOUT)
  }

  runtime {
    preemptible:      10
    docker:           docker_image
    memory:           memory
    disks:            storage
    ignoreExitStatus: true
  }
}


task  CollectShards {
  input {
    String          rundir              # url to a run directory
    Array[Object]   non_compliant       # each object in this array has
                                        # members path and type
    Array[Object]   sequencing_sentinel # ensures proper sequencing of
                                        # tasks; it is otherwise ignored
    String          docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"
  String   DUMMY     = OUTPUTDIR + "/DUMMY"
  String   DEBUGGING = OUTPUTDIR + "/DEBUGGING"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  mkdir --parents '~{OUTPUTDIR}'
  ( hostname --long ; date ) > '~{DUMMY}'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  printf -- '{}' > '~{DEBUGGING}'

  export COLLECT_SHARDS__DEBUGGING='~{DEBUGGING}'

  collect_shards.py                  \
      --logging-level=DEBUG          \
      '~{write_json(non_compliant)}' \
      '~{rundir}'                    \
    | tee '~{STDOUT}'

  date >> '~{DUMMY}'
  >>>

  output {
    String                 sequencing_dummy = read_string(DUMMY)
    Array[Object]   storage_and_batch_array = read_json(STDOUT)
    File                     debugging_info = DEBUGGING
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}


task  ConcatenateShards {
  # ...
  input {
    Int             storage        # required disk size (in GB)
    String          rundir         # url to a run directory
    Array[Object]   batch          # each object in this array has
                                   # members nshards, source, and
                                   # endpoint
    String          docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"
  String   DUMMY     = OUTPUTDIR + "/DUMMY"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  mkdir --parents '~{OUTPUTDIR}'
  ( hostname --long ; date ) > '~{DUMMY}'

  concatenate_shards.py      \
      --logging-level=DEBUG  \
      '~{rundir}'            \
      '~{write_json(batch)}' \
    | tee '~{STDOUT}'

  date >> '~{DUMMY}'
  >>>

  output {
    # String   sequencing_dummy = read_string(DUMMY)
    Object   results = read_json(STDOUT)
  }

  runtime {
    preemptible:      3
    docker:           docker_image
    disks:            "local-disk ~{storage} HDD"
    ignoreExitStatus: true
  }
}


task CheckBatchedResults {
  input {
    Array[Object]   batched_results
    String          docker_image
  }

  String   OUTPUTDIR  = "OUTPUT"
  String   STDOUT     = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  RESULTSDIR="${PWD}/results"

  mkdir --parents "${RESULTSDIR}"

  RESULTS="${RESULTSDIR}/results.json"
  consolidate_batched_results.py       \
      '~{write_json(batched_results)}' \
    > "${RESULTS}"

  mkdir --parents '~{OUTPUTDIR}'

  if results_ok.py "${RESULTS}"
  then
      ISOK=true
  else
      ISOK=false
  fi

  printf -- '%s\n' "${ISOK}" | tee '~{STDOUT}'

  >>>

  output {
    Boolean isok = read_json(STDOUT)
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}


task  MakePushBatches {

  input {
    String   rundir
    String   docker_image
  }

  String   OUTPUTDIR  = "OUTPUT"
  String   STDOUT     = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  mkdir --parents '~{OUTPUTDIR}'

  make_push_batches.py '~{rundir}' | tee '~{STDOUT}'

  >>>

  output {
    Array[Array[String]]   batches = read_json(STDOUT)
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}

task  PushScrubbed {

  input {
    String          rundir
    String          release_datadir
    Array[String]   batch
    String          docker_image
  }

  String   OUTPUTDIR  = "OUTPUT"
  String   STDOUT     = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{rundir}'

  /mgbpmbiofx/packages/biofx-orchestration-utils/bin/setup-rclone-remote.sh -p mgb-lmm-gcp-infrast-1651079146 -w prod-biobank-scrub -r '~{release_datadir}'

  mkdir --parents '~{OUTPUTDIR}'

  push_scrubbed.py           \
      '~{rundir}'            \
      '~{release_datadir}'   \
      '~{write_json(batch)}' \
    | tee '~{STDOUT}'

  >>>

  output {
    Object   results = read_json(STDOUT)
  }

  runtime {
    preemptible:      10
    docker:           docker_image
    ignoreExitStatus: true
  }
}

task  Summarize {

  input {
    Array[Object]   scrub_results
    Array[Object]   concatenation_results
    Array[Object]   push_results
    String          docker_image
  }

  String   OUTPUTDIR  = "OUTPUT"
  String   STDOUT     = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  RESULTSDIR="${PWD}/results"

  mkdir --parents "${RESULTSDIR}"

  SCRUBRESULTS="${RESULTSDIR}/scrub.json"
  consolidate_batched_results.py     \
      '~{write_json(scrub_results)}' \
    > "${SCRUBRESULTS}"

  CONCATENATIONRESULTS="${RESULTSDIR}/concatenation.json"
  consolidate_batched_results.py             \
      '~{write_json(concatenation_results)}' \
    > "${CONCATENATIONRESULTS}"

  PUSHRESULTS="${RESULTSDIR}/push.json"
  consolidate_batched_results.py     \
      '~{write_json(push_results)}' \
    > "${PUSHRESULTS}"

  mkdir --parents '~{OUTPUTDIR}'

  summarize.py                  \
      "${SCRUBRESULTS}"         \
      "${CONCATENATIONRESULTS}" \
      "${PUSHRESULTS}"          \
    | tee '~{STDOUT}'

  >>>

  output {
    File   results = STDOUT
  }

  runtime {
    preemptible: 2
    docker:      docker_image
  }
}

# -----------------------------------------------------------------------------
