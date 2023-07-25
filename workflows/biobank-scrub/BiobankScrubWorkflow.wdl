version 1.0


workflow BiobankScrubWorkflow {
  input {
    String         docker_image         = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/biobank_scrubbing:1.0.0"

    Int            nbatches             = 500
    String         initial_datadir      = "gs://url/to/default/initial/datadir"
    String         current_datadir
    String         staging_area
    File?          maybe_withdrawn_list
  }

  if (!defined(maybe_withdrawn_list)) {
    call FetchWithdrawnList {
      input:
        docker_image = docker_image
    }
  }

  File withdrawn_list = select_first([
    maybe_withdrawn_list,
    FetchWithdrawnList.withdrawn_list
  ])

  call GetDatasetIds {
    input:
      base_datadir = current_datadir,
      docker_image = docker_image
  }

  scatter (dataset_id in GetDatasetIds.dataset_ids) {
    call FindNonCompliant {
      input:
        dataset_id      = dataset_id,
        withdrawn_list  = withdrawn_list,
        current_datadir = current_datadir,
        initial_datadir = initial_datadir,
        docker_image    = docker_image
    }
  }

  Array[Object] non_compliant = flatten(FindNonCompliant.non_compliant)

  call MakeBatches {
    input:
      non_compliant = non_compliant,
      nbatches      = nbatches,
      staging_area  = staging_area,   # to check for already-scrubbed files
      docker_image  = docker_image
  }

  scatter (batch in MakeBatches.batches) {
    call ScrubBatch {
      input:
        withdrawn_list  = withdrawn_list,
        initial_datadir = initial_datadir,
        staging_area    = staging_area,
        batch           = batch,
        docker_image    = docker_image
    }
  }

  # -------------------------------------------------------------------

  call CollectShards {
    input:
      staging_area    = staging_area,
      non_compliant   = non_compliant,
      scrub_results   = flatten(ScrubBatch.results),
      docker_image    = docker_image
  }

  scatter (shards in CollectShards.shards) {
    call ConcatenateShards {
      input:
        staging_area = staging_area,
        shards       = shards,
        docker_image = docker_image
    }
  }

  call VerifyRelease {
    input:
      sequencing_dummy = ScrubBatch.sequencing_dummy,
      release_dir      = staging_area,
      non_compliant    = non_compliant,
      docker_image     = docker_image
  }

  output {
    File   scrub_results       = write_json(flatten(ScrubBatch.results))
    File   verification_report = write_json(VerifyRelease.report)
  }
}

# -----------------------------------------------------------------------------

task  FetchWithdrawnList {
  input {
    String  docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  mkdir --parents '~OUTPUTDIR'
  list_withdrawn_subjects.py > '~STDOUT'
  >>>

  output {
    File   withdrawn_list = STDOUT
  }

  runtime {
    docker: docker_image
  }
}


task  GetDatasetIds {
  input {
    String  base_datadir
    String  docker_image
  }

  String   OUTPUTDIR = "OUTPUT"
  String   STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  # set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  # list() {
  #     local directory="${1}"
  #     if [[ ${directory} == 'gs://'* ]]
  #     then
  #         gsutil -m ls "${directory}"
  #     else
  #         find "${directory}" -mindepth 1 -maxdepth 1
  #     fi
  # }
  #
  # mkdir --parents '~OUTPUTDIR'
  # list '~{base_datadir}' | grep -P '/\d{4}/$' | sort > '~STDOUT'

  mkdir --parents '~OUTPUTDIR'
  list_dataset_ids.py '~{base_datadir}' > '~STDOUT'
  >>>

  output {
    Array[String]   dataset_ids = read_lines(STDOUT)
  }

  runtime {
    docker: docker_image
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

  mkdir --parents '~OUTPUTDIR'
  find_non_compliant.py    \
      '~{withdrawn_list}'  \
      '~{current_datadir}' \
      '~{initial_datadir}' \
      '~{dataset_id}'      \
    > '~{STDOUT}'
  >>>

  output {
    Array[Object]   non_compliant = read_json(STDOUT)
  }

  runtime {
    docker: docker_image
  }
}

task  MakeBatches {
  input {
    Array[Object]   non_compliant   # each object in this array has
                                    # members path and type
    Int             nbatches        # desired number of batches
    String          staging_area    # url to a directory for staging
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

  mkdir --parents '~OUTPUTDIR'
  make_batches.py                    \
      '~{write_json(non_compliant)}' \
      ~{nbatches}                    \
      '~{staging_area}'              \
    > '~{STDOUT}'
  # PRINTS TO STDOUT AN ARRAY OF nbatches ARRAYS OF OBJECTS; EACH OBJECT IN
  # THIS OUTPUT CORRESPONDS TO A (STILL-PENDING) OBJECT IN THE non_compliant
  # INPUT.
  >>>

  output {
    Array[Array[Object]]   batches = read_json(STDOUT)
  }

  runtime {
    docker: docker_image
  }
}

task  ScrubBatch {
  input {
    File            withdrawn_list
    String          initial_datadir
    String          staging_area
    Array[Object]   batch
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

  ( hostname --long ; date ) > '~{DUMMY}'

  scrub_batch.py             \
      '~{withdrawn_list}'    \
      '~{initial_datadir}'   \
      '~{staging_area}'      \
      '~{write_json(batch)}' \
    > '~{STDOUT}'

  date >> '~{DUMMY}'
  >>>

  output {
    String          sequencing_dummy = read_lines(DUMMY)
    Array[Object]   results          = read_json(STDOUT)
  }

  runtime {
    docker: docker_image
    memory: "100GB"
    cpus:   8
    disks:  "local-disk 750 HDD" ### FIXME: find smallest viable disk size
  }
}

task CollectShards {
  input {
    String          staging_area
    Array[Object]   non_compliant
    Array[Object]   scrub_results
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

  ( hostname --long ; date ) > '~{DUMMY}'

  collect_shards.py                  \
      '~{staging_area}'              \
      '~{write_json(non_compliant)}' \
      '~{write_json(scrub_results)}' \
    > '~{STDOUT}'

  date >> '~{DUMMY}'
  >>>

  output {
    String                 sequencing_dummy = read_lines(DUMMY)
    Array[Array[Object]]   shards           = read_json(STDOUT)
  }

  runtime {
    docker: docker_image
    memory: "100GB"
    cpus:   8
    disks:  "local-disk 750 HDD" ### FIXME: find smallest viable disk size
  }
}

task  VerifyRelease {
  input {
    Array[String]   sequencing_dummy   # ignored; needed only to ensure proper
                                       # sequencing of tasks
    String          release_dir        # url to a release's root
    Array[Object]   non_compliant      # array of non-compliant items
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

  mkdir --parents '~OUTPUTDIR'
  verify_release.py '~{release_dir}' '~{write_json(non_compliant)}' > '~{STDOUT}'
  >>>

  output {
    File   report = STDOUT
  }

  runtime {
    docker: docker_image
  }
}
