version 1.0
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# -----------------------------------------------------------------------------

workflow SubsetShardsWorkflow {
  input {
    String   sourcebase
    String?  shardsbase
    String   targetbase
    File     shardsmapfile
    String   workspace
    Boolean  noconcat      = false
    Boolean? nodelete
    File?    regions
    File?    targets
    File?    samples
    File?    script
    Int      nbatches      = 500
    String   docker_image  = "us-central1-docker.pkg.dev/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:latest"
  }

  if (defined(shardsbase)) {
    call FetchSentinels {
      input:
        docker_image = docker_image
      , shardsbase   = select_first([shardsbase])
      , workspace    = workspace
    }
  }

  call CheckInputsResolveShardsbaseBatchShards {
    input:
      docker_image  = docker_image
    , noregions     = !defined(regions)
    , notargets     = !defined(targets)
    , nosamples     = !defined(samples)
    , sourcebase    = sourcebase
    , shardsbase    = shardsbase
    , targetbase    = targetbase
    , noconcat      = noconcat
    , nodelete      = nodelete
    , sentinels     = select_first([FetchSentinels.sentinels, []])
    , shardsmapfile = shardsmapfile
    , nbatches      = nbatches
  }

  Array[Array[String]] batches = CheckInputsResolveShardsbaseBatchShards.batches
  Object               params  = CheckInputsResolveShardsbaseBatchShards.params

  scatter (batch in batches) {
    call SubsetShards {
      input:
        docker_image = docker_image
      , script       = script
      , regions      = regions
      , targets      = targets
      , samples      = samples
      , sourcebase   = sourcebase
      , shardsbase   = params.shardsbase
      , batch        = batch
      , workspace    = workspace
    }
  }

  Boolean doconcat = params.doconcat

  if (doconcat) {

    Array[Object] shardsmap = CheckInputsResolveShardsbaseBatchShards.shardsmap
    String        prefix    = params.shardsbase

    scatter (spec in shardsmap) {

      String vcf     = spec.vcf
      String basedir = "~{prefix}/data/~{vcf}" # sic!

      call GetTotalSize {
        input:
          docker_image = docker_image
        , basedir      = basedir
        , workspace    = workspace
        , sequencing   = SubsetShards.sequencing
      }

      call ConcatenateShards {
        input:
          docker_image = docker_image
        , storage      = 3 * GetTotalSize.footprint + 10
        , basedir      = basedir
        , shards       = spec.shards
        , target       = "~{targetbase}/~{vcf}"
        , workspace    = workspace
      }
    }

    Boolean dodelete = params.dodelete

    if (dodelete) {
      call DeleteShards {
        input:
          docker_image = docker_image
        , shardsbase   = params.shardsbase
        , workspace    = workspace
        , sequencing   = ConcatenateShards.sequencing
      }
    }
  }

}


task FetchSentinels {
  input {
    String docker_image
    String shardsbase
    String workspace
  }

  String SENTINELS = "output/SENTINELS"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # export RCLONE_LOG_LEVEL=DEBUG
  # export RCLONE_STATS_LOG_LEVEL=DEBUG

  typeset -p >&2
  rclone version >&2

  mkdir --verbose --parents "$( dirname '~{SENTINELS}' )"

  export WORKSPACE='~{workspace}'
  SHARDSBASE="$( mapurl.sh '~{shardsbase}' )"

  rclone                        \
      lsf                       \
      --recursive               \
      --files-only              \
      "${SHARDSBASE}/sentinels" \
    > '~{SENTINELS}'

  >>>

  output {
    Array[String] sentinels = read_lines(SENTINELS)
  }

  runtime {
    preemptible: 5
    docker     : docker_image
  }
}

task CheckInputsResolveShardsbaseBatchShards {
  input {
    String        docker_image
    Boolean       noregions
    Boolean       notargets
    Boolean       nosamples
    String        sourcebase
    String?       shardsbase
    String        targetbase
    Boolean       noconcat
    Boolean?      nodelete
    Array[String] sentinels
    File          shardsmapfile
    Int           nbatches
    # -----------------------------------------------------------------
    Boolean       noshuffle     = false
    Int?          seed
  }

  Boolean no_sentinels = length(sentinels) == 0

  String  BATCHES      = "output/BATCHES"
  String  PARAMS       = "output/PARAMS"

  command <<<
  free --bytes >&2
  free --human >&2

  mkdir --verbose --parents "$( dirname '~{BATCHES}'   )"
  mkdir --verbose --parents "$( dirname '~{PARAMS}'    )"

  python3 <<EOF
  import sys
  import os
  import random
  import json
  import uuid
  import re

  def error(message):
      print(f'ERROR: {message}', file=sys.stderr)
      sys.exit(1)


  def unsafe_endpoint(endpoint, sourcebase='~{sourcebase}'.rstrip('/')):
      return f'{endpoint.rstrip("/")}/'.startswith(f'{sourcebase}/')


  def checkinputs():

      missing_inputs = ~{if noregions && notargets && nosamples
                         then "True" else "False"}
      if missing_inputs:
          error('Neiter regions nor targets nor samples are specified')

      targetbase = '~{targetbase}'.rstrip('/')
      if unsafe_endpoint(targetbase):
          error('writing to targetbase would modify sourcebase')

      noconcat = ~{if noconcat then "True" else "False"}

      shardsbase_defined = ~{if defined(shardsbase) then "True" else "False"}

      if shardsbase_defined and unsafe_endpoint('~{shardsbase}'):
          error('writing to shardsbase would modify sourcebase')

      if noconcat:
          if ~{if defined(nodelete) && !select_first([nodelete])
               then "True" else "False"}:
              error('inconsistent inputs: noconcat=true and nodelete=false')
          if shardsbase_defined:
              error('inconsistent inputs: noconcat=true and '
                    'shardsbase specified')
          shardsbase = targetbase
          nodelete = True
      else:
          nodelete = ~{if select_first([nodelete, false])
                       then "True" else "False"}

          if shardsbase_defined:
              shardsbase = '~{shardsbase}'
          else:
              tmpdir = os.path.join(re.match(r'^[^:]+:/*[^/]+',
                                    '~{targetbase}').group(0),
                                    'tmp')
              # FIXME: the use of uuid.uuid4 here breaks call-caching
              shardsbase = os.path.join(tmpdir, str(uuid.uuid4()))

              if nodelete:
                  print( 'INFO: shards will be kept under '
                        f'\047{shardsbase}\047',
                        file=sys.stderr)

      params = dict(doconcat=not noconcat,
                    dodelete=not nodelete,
                    shardsbase=shardsbase)

      with open('~{PARAMS}', 'w') as writer:
          json.dump(params, writer)


  def readjson(filepath):
      with open(filepath) as reader:
          return json.load(reader)


  def batch(cases, nbatches):

      ncases = len(cases)
      min_stride, leftover = divmod(ncases, nbatches)
      batches = []
      offset = 0
      while offset < ncases:
          stride = min_stride
          if leftover > 0:
              stride += 1
              leftover -= 1
          elif stride == 0:
              assert ncases < nbatches
              break

          batches.append(cases[offset:offset + stride])
          offset += stride

      return batches


  def check_shardsmap(shardsmap):
      if len(shardsmap) == 0:
          error('empty shardsmap')

      seen = set()
      for spec in shardsmap:
          vcf = spec['vcf']
          if vcf == '':
              error('invalid shardsmap: empty vcf value')
          if vcf in seen:
              error('invalid shardsmap: repeated vcf value')
          seen.add(vcf)

          prefix = f'invalid shardsmap for vcf={vcf!r}'
          shards = spec['shards']
          nshards = len(shards)
          if nshards == 0:
              error(f'{prefix}: empty shards list')
          elif '' in shards:
              error(f'{prefix}: shards list contains empty paths')
          elif nshards > len(set(shards)):
              error(f'{prefix}: shards list contains repeats')


  def build_cases(shardsmap):

      check_shardsmap(shardsmap)

      ~{if no_sentinels
        then "# NB: the sentinels variable will not be used in this execution"
        else ""}
      if ~{if no_sentinels then "True" else "False"}:
          sentinels = None
      else:
          sentinels = set(['~{sep="', '" sentinels}'])

      cases = [f"{relpath}\t{spec['vcf']}"
                   for spec in shardsmap
                       for relpath in spec['shards']
               ~{if no_sentinels then ""
                 else "if not relpath in sentinels"}]

      if ~{if noshuffle then "False" else "True"}:
          ~{"random.seed(" + seed + ")"}
          random.shuffle(cases)

      return cases


  def main():

      checkinputs()

      shardsmap = readjson('~{shardsmapfile}')

      cases = build_cases(shardsmap)

      batches = batch(cases, ~{nbatches})

      with open('~{BATCHES}', 'w') as writer:
          json.dump(batches, writer)


  # --------------------------------------------------------------------------

  main()
  EOF
  >>>

  output {
    Array[Object]        shardsmap = read_json(shardsmapfile)
    Array[Array[String]] batches   = read_json(BATCHES)
    Object               params    = read_json(PARAMS)
  }

  runtime {
    preemptible: 5
    docker     : docker_image
  }
}


task SubsetShards {
  input {
    String        docker_image
    File?         script
    File?         regions
    File?         targets
    File?         samples
    String        sourcebase    
    String        shardsbase
    Array[String] batch
    String        workspace
    # ---------------------------------------------------------------
    Int           preemption   = 10
  }

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  export RCLONE_LOG_LEVEL=DEBUG
  export RCLONE_STATS_LOG_LEVEL=DEBUG

  typeset -p >&2
  compgen -abckA function | sort --unique >&2
  rclone version >&2

  fileexists() {
      local url="${1%/}"
      local parent="$( dirname "${url}" )"
      local basename="$( basename "${url}" )"
      local -a found
      readarray -t found < <(
                              rclone                      \
                                  lsf                     \
                                  --files-only            \
                                  --include="${basename}" \
                                  "${parent}"
                            )

      local nfound="${#found[@]}"
      if (( nfound == 0 ))
      then
          printf -- 'false'
      elif (( nfound == 1 ))
      then
          printf -- 'true'
      else
          printf -- 'fileexists: unexpected nfound=%d\n' "${nfound}" >&2
          return 1
      fi
  }

  footprint() {
      (
        set +o xtrace
        set +o errexit
        set +o pipefail
        set +o nounset

        exec >&2

        printf -- '\n\n'
        df --block-size=1   / /* ${TMPDIR} 2>/dev/null
        printf -- '\n'
        df --human-readable / /* ${TMPDIR} 2>/dev/null
        printf -- '\n\n'

        du --summarize --one-file-system --apparent-size --block-size=1   / /* ${TMPDIR} 2>/dev/null
        printf -- '\n'
        du --summarize --one-file-system --apparent-size --human-readable / /* ${TMPDIR} 2>/dev/null
        printf -- '\n\n'
      )
  }

  export WORKSPACE='~{workspace}'
  SOURCEBASE="$( mapurl.sh '~{sourcebase}' )"
  SHARDSBASE="$( mapurl.sh '~{shardsbase}' )"

  SUBSETTING_FLAGS=()

  if ~{if defined(regions) then "true" else "false"}
  then
      SUBSETTING_FLAGS+=( --regions-file '~{regions}' )
  fi

  if ~{if defined(targets) then "true" else "false"}
  then
      SUBSETTING_FLAGS+=( --targets-file '~{targets}' )
  fi

  if ~{if defined(samples) then "true" else "false"}
  then
      SUBSETTING_FLAGS+=( --samples-file '~{samples}' )
  fi

  if '~{if defined(script) then "true" else "false"}'
  then
      SCRIPT='~{script}'
      chmod +x "${SCRIPT}"
  else
      SCRIPT=subsetvcf.sh
  fi

  BATCH='~{write_lines(batch)}'

  footprint

  while IFS=$'\t' read -d$'\n' RELPATH SUBDIR
  do
      SENTINEL="${SHARDSBASE}/sentinels/${RELPATH}"
      if "$( fileexists "${SENTINEL}" )"
      then
          (
            set +o xtrace
            exec >&2
            printf -- '\nSubsetShards: INFO: '
            printf -- 'skipping \047%s\047 (already done)\n\n' \
                      "${RELPATH}"
          )
          continue
      fi

      INPUTSHARD="$( mktemp )"
      OUTPUTSHARD="$( mktemp )"

      rclone copyto "${SOURCEBASE}/${RELPATH}"     "${INPUTSHARD}"
      rclone copyto "${SOURCEBASE}/${RELPATH}.tbi" "${INPUTSHARD}.tbi"

      touch --date '1 Jan 1970 00:00:00 -0000' "${INPUTSHARD}"
      touch --date '1 Jan 1970 00:00:01 -0000' "${INPUTSHARD}.tbi"

      "${SCRIPT}"                                         \
          "${INPUTSHARD}"                                 \
          "${OUTPUTSHARD}"                                \
          ${SUBSETTING_FLAGS[@]+"${SUBSETTING_FLAGS[@]}"}

      SHARDPATH="${SUBDIR}/${RELPATH}"
      rclone copyto "${OUTPUTSHARD}"     "${SHARDSBASE}/data/${SHARDPATH}"
      rclone copyto "${OUTPUTSHARD}.tbi" "${SHARDSBASE}/data/${SHARDPATH}.tbi"

      rclone touch "${SENTINEL}"

      rm "${INPUTSHARD}"  "${INPUTSHARD}.tbi"  \
         "${OUTPUTSHARD}" "${OUTPUTSHARD}.tbi"

      footprint

  done < "${BATCH}"
  >>>

  output {
    Boolean sequencing = true
  }

  runtime {
    preemptible: preemption
    docker     : docker_image
  }
}

task GetTotalSize {
  input {
    String         docker_image
    String         basedir
    String         workspace
    Array[Boolean] sequencing
  }

  String SIZE = "output/SIZE"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # export RCLONE_LOG_LEVEL=DEBUG
  # export RCLONE_STATS_LOG_LEVEL=DEBUG

  typeset -p >&2
  rclone version >&2

  export WORKSPACE='~{workspace}'
  BASEDIR="$( mapurl.sh '~{basedir}' )"

  mkdir --verbose --parents "$( dirname '~{SIZE}' )"

  BYTES="$(
            rclone size "${BASEDIR}"               \
              | tail --lines=1                     \
              | perl -lpe 's/^.*?(\d+) Byte.*/$1/'
          )"

  perl -e "printf qq(%f\n), ${BYTES}/2**30" > '~{SIZE}'
  printf -- 'Size found for \047%s\047: %.1f GiB\n' '~{basedir}' "$( cat '~{SIZE}' )"
  >>>

  output {
    Int footprint = ceil(read_float(SIZE))
  }

  runtime {
    preemptible: 5
    docker     : docker_image
  }
}


task ConcatenateShards {
  input {
    String        docker_image
    Int           storage
    String        basedir
    Array[String] shards
    String        target
    String        workspace

    Int           memory       = 16
  }

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  export RCLONE_LOG_LEVEL=DEBUG
  export RCLONE_STATS_LOG_LEVEL=DEBUG

  typeset -p >&2
  rclone version >&2

  free --bytes >&2
  free --human >&2

  export WORKSPACE='~{workspace}'
  BASEDIR="$( mapurl.sh '~{basedir}' )"
  TARGET="$( mapurl.sh '~{target}' )"

  LOCALBASEDIR="$( mktemp --directory )"

  rclone copy "${BASEDIR}" "${LOCALBASEDIR}"

  find "${LOCALBASEDIR}" -type f

  SHARDS="$( mktemp )"
  perl -lpe "s@^@${LOCALBASEDIR}/@" '~{write_lines(shards)}' > "${SHARDS}"

  LOCALTARGET="$( mktemp )"

  # --------------------------------------------------------------------------

  /usr/bin/time --verbose         \
      bcftools                    \
          concat                  \
          --file-list="${SHARDS}" \
          --no-version            \
          --naive                 \
        > "${LOCALTARGET}"

  /usr/bin/time --verbose  \
      bcftools             \
          index            \
          --force          \
          --tbi            \
          "${LOCALTARGET}"

  # --------------------------------------------------------------------------

  for EXTENSION in '' '.tbi'
  do
      rclone copyto "${LOCALTARGET}${EXTENSION}" "${TARGET}${EXTENSION}"
  done
  >>>

  output {
    Boolean sequencing = true
  }

  runtime {
    preemptible: 3
    docker     : docker_image
    disks      : "local-disk ~{storage} HDD"
    memory     : "~{memory} GB"
  }
}


task DeleteShards {
  input {
    String         docker_image
    String         shardsbase
    String         workspace
    Array[Boolean] sequencing
  }

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # export RCLONE_LOG_LEVEL=DEBUG
  # export RCLONE_STATS_LOG_LEVEL=DEBUG

  typeset -p >&2
  rclone version >&2

  export WORKSPACE='~{workspace}'
  SHARDSBASE="$( mapurl.sh '~{shardsbase}' )"

  rclone purge "${SHARDSBASE}"
  >>>

  runtime {
    preemptible: 5
    docker     : docker_image
  }
}
