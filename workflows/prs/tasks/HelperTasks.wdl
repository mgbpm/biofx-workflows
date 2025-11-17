version 1.0

task GetBaseMemory {
  input {
    File?  vcf
    Int?   nvariants
    String docker_image = "python:3.11"
    Int    addldisk     = 20
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  Int     file_size         = 2 * ceil(size(vcf, "GB"))
  Int     final_disk_size   = addldisk + file_size
  Boolean ERROR             = defined(vcf) == defined(nvariants)
  String  OUTPUTDIR         = "OUTPUT"
  String  NVARIANTS         = OUTPUTDIR + "/nvariants.txt"
  String  GIGABYTES         = OUTPUTDIR + "/gigabytes.txt"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  if ~{if ERROR then "true" else "false"}
  then
      printf -- 'INTERNAL ERROR: too few or too many arguments specified' >&2
      exit 1
  fi

  # --------------------------------------------------------------------------

  mkdir --verbose --parents '~{OUTPUTDIR}'

  NVARIANTS=~{if defined(nvariants)
              then nvariants
              else "\"$( zgrep --count --invert-match '^#' '" + vcf + "' | tee '" + NVARIANTS + "' )\""}

  python3 <<EOF > '~{GIGABYTES}'
  import math
  print(8 + max(0, math.ceil((${NVARIANTS} - 50000000)/10000000)))
  EOF
  >>>

  output {
    Int gigabytes  = read_int(GIGABYTES)
    Int nvariants_ = if defined(nvariants) then nvariants else read_int(NVARIANTS)
  }

  runtime {
    docker: "~{docker_image}"
    disks: "local-disk ~{final_disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task RenameChromosomesInTsv {
  input {
    File    tsv
    Boolean skipheader
    File    lookup      = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
    String docker_image = "python:3.11"
    Int    addldisk     = 20
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  Int    file_size       = 2 * ceil(size(tsv, "GB"))
  Int    final_disk_size = addldisk + file_size
  String RENAMED         = "OUTPUT/renamed/" + basename(tsv)

  command <<<
  python3 <<EOF
  import sys
  import os
  import re

  def error(message):
      print(f'ERROR: {message}', file=sys.stderr)
      sys.exit(1)


  def read_lookup():
      lookup = dict()
      with open('~{lookup}') as reader:
          for rawline in reader:
              key, value = rawline.rstrip('\r\n').split('\t')
              lookup[key] = value
      return lookup


  def main():

      def rename(chromosomename,
                 _lookup=read_lookup(),
                 _parse_re=re.compile(r'^([\da-z]+)(.*)',
                                      flags=re.I)):

          match = _parse_re.search(chromosomename)

          if match:
              core, rest = match.groups()
              if core in _lookup:
                  return f'{_lookup[core]}{rest}'

          error(f'Unsupported H. sapiens chromosome name: {chromosomename}')

      inputtsv = '~{tsv}'
      outputtsv = '~{RENAMED}'

      os.makedirs(os.path.dirname(outputtsv), exist_ok=True)

      with open(outputtsv, 'w') as writer:

          with open(inputtsv) as reader:

              skipheader = ~{if skipheader then "True" else "False"}

              if skipheader:
                  writer.write(next(reader))

              for rawline in reader:
                  row = rawline.rstrip('\r\n').split('\t')
                  parts = row[0].split(':')
                  newid = ':'.join([rename(parts[0])] + parts[1:])
                  print('\t'.join([newid] + row[1:]), file=writer)
                  continuing = True

  # --------------------------------------------------------------------------

  main()
  EOF
  >>>

  output {
    File renamed = RENAMED
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{final_disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task RenameChromosomesInVcf {
  input {
    File   vcf
    File   lookup       = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
    String docker_image = "biocontainers/bcftools:v1.9-1-deb_cv1"
    Int    addldisk     = 20
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  Int    file_size       = 2 * ceil(size(vcf, "GB"))
  Int    final_disk_size = addldisk + file_size
  String RENAMED = "OUTPUT/renamed/" + basename(vcf)

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents "$( dirname '~{RENAMED}' )"

  WORKDIR="$( mktemp --directory )"
  INPUTVCF="${WORKDIR}/input.vcf.gz"

  ln --symbolic --verbose '~{vcf}' "${INPUTVCF}"

  bcftools index    \
      --force       \
      --tbi         \
      "${INPUTVCF}"

  bcftools annotate                    \
      --no-version                     \
      --output='~{RENAMED}'            \
      --output-type=z                  \
      --rename-chr='~{lookup}'         \
      --set-id '%CHROM:%POS:%REF:%ALT' \
      "${INPUTVCF}"

  >>>

  output {
    File renamed = RENAMED
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{final_disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task SubsetVcf {
  input {
    File    inputvcf
    File    regions
    String  label        = "data"
    Boolean nocleanup    = false
    String  docker_image = "biocontainers/bcftools:v1.9-1-deb_cv1"
    Int     addldisk     = 20
    Int     mem_size     = 2
    Int     preemptible  = 1
  }

  Int    file_size       = 3 * ceil(size(inputvcf, "GB"))
  Int    final_disk_size = addldisk + file_size
  String OUTPUTDIR       = "OUTPUT"
  String OUTPUTVCF       = OUTPUTDIR + "/" + label + ".vcf.gz"
  String NREGIONS        = OUTPUTDIR + "/NREGIONS"

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{final_disk_size}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents '~{OUTPUTDIR}'
  WORKDIR="$( mktemp --directory )"
  INPUTVCF="${WORKDIR}/input.vcf.gz"

  ln --symbolic --verbose '~{inputvcf}' "${INPUTVCF}"

  bcftools index      \
      --force         \
      --tbi           \
      "${INPUTVCF}"

  cleanup() {

      bcftools                             \
          norm                             \
          --multiallelics -any             \
          --no-version                     \
          --output-type v                  \
    | bcftools                             \
          annotate                         \
          --no-version                     \
          --output-type v                  \
          --remove 'INFO,FORMAT'           \
          --set-id '%CHROM:%POS:%REF:%ALT'

  }

  if ~{if nocleanup then "true" else "false"}
  then
      POSTPROCESS=cat
  else
      POSTPROCESS=cleanup
  fi

  bcftools                             \
      view                             \
      --no-version                     \
      --output-type v                  \
      --regions-file '~{regions}'      \
      "${INPUTVCF}"                    \
    | "${POSTPROCESS}"                 \
    | bcftools                         \
          view                         \
          --no-version                 \
          --output-type z              \
          --output-file '~{OUTPUTVCF}'

  wc --lines < '~{regions}' > '~{NREGIONS}'

  # ---------------------------------------------------------------------------

  printf -- 'FINAL STORAGE UTILIZATION:\n'
  df --human

  # ---------------------------------------------------------------------------
  >>>

  output {
    File result   = OUTPUTVCF
    Int  nregions = read_int(NREGIONS)
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{final_disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task Union {
  input {
    Array[File]+ lists
    String       docker_image = "ubuntu:21.10"
    Int          addldisk     = 20
    Int          mem_size     = 2
    Int          preemptible  = 1
  }
  
  Int    final_disk_size = 20 + addldisk
  String OUTPUTDIR       = "OUTPUT"
  String RESULT          = OUTPUTDIR + "/RESULT"

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset
    # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
    # set -o xtrace

    mkdir --verbose --parents '~{OUTPUTDIR}'

    sort --unique '~{sep="' '" lists}' > '~{RESULT}'
  >>>

  output {
    File result = RESULT
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{final_disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task Intersection {
  input {
    Array[File]+ lists
    String       docker_image = "ubuntu:21.10"
    Int          addldisk     = 20
    Int          mem_size     = 2
    Int          preemptible  = 1
  }

  Int    final_disk_size = 20 + addldisk
  String OUTPUTDIR       = "OUTPUT"
  String RESULT          = OUTPUTDIR + "/RESULT"

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset
    # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
    # set -o xtrace

    INPUTS=( '~{sep="' '" lists}' )

    HOLD="$( mktemp )"

    mkdir --verbose --parents '~{OUTPUTDIR}'

    sort --unique "${INPUTS[0]}" > '~{RESULT}'
    for INPUT in "${INPUTS[@]:1}"
    do
        comm -1 -2 '~{RESULT}' <( sort --unique "${INPUT}" ) > "${HOLD}"
        mv "${HOLD}" '~{RESULT}'
    done
  >>>

  output {
    File result = RESULT
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{final_disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task ListShards {
  input {
    String source
    String workspace
    String docker_image
    Int    disk_size    = 10
    Int    mem_size     = 2
    Int    preemptible  = 1
  }

  String OUTPUTDIR = "OUTPUT"
  String STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset
    set -o xtrace
    # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  
    # -----------------------------------------------------------------------------
    export WORKSPACE='~{workspace}'
    SOURCE="$( mapurl.sh '~{source}' )"
  
    mkdir --parents '~{OUTPUTDIR}'
  
    rclone lsf --files-only --recursive "${SOURCE}" \
      | grep --perl-regexp 'shards/.*\.gz$'         \
      | perl -lne '
          BEGIN {
            %lookup = (
                        X  => 23,
                        Y  => 24,
                        XY => 25,
                        MT => 26
                      )
          }
          $_ =~ /chr(\d+|XY|X|Y|MT)\b/;
          $index = ( $lookup{$1} or $1 );
          printf qq(%s\t%s\n), $index, $_;
         '                                          \
      | sort                                        \
            --field-separator=$'\t'                 \
            --key=1,1n                              \
            --key=2,2V                              \
      | cut --fields=2                              \
      | tee '~{STDOUT}'                             \
      >&2
  >>>

  output {
    File relpaths = STDOUT
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task MakeBatches {
  input {
    File          cases
    Int           nbatches
    Array[String] exclude   = []
    Boolean       noshuffle = false
    Int?          seed
    String        docker_image = "python:3.9.10"
    Int           disk_size    = 10
    Int           mem_size     = 2
    Int           preemptible  = 1
  }

  Boolean no_exclude = length(exclude) == 0

  String OUTPUTDIR     = "OUTPUT"
  String BATCHES       = OUTPUTDIR + "/BATCHES"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  mkdir --parents '~{OUTPUTDIR}'

  python3 <<EOF
  import sys
  import os
  import random
  import json


  def getcases():

      ~{if no_exclude
        then "# NB: the exclude variable will not be used in this execution"
        else ""}
      if ~{if no_exclude then "True" else "False"}:
          exclude = None
      else:
          exclude = set(['~{sep="', '" exclude}'])

      return [case for case in readlines('~{cases}')
              ~{if no_exclude then "" else "if not case in exclude"}]


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


  def readlines(filepath):
      with open(filepath) as reader:
          return [rawline.rstrip('\n') for rawline in reader]


  def main():

      cases = getcases()

      if ~{if noshuffle then "False" else "True"}:
          ~{"random.seed(" + seed + ")"}
          random.shuffle(cases)

      batches = batch(cases, ~{nbatches})

      with open('~{BATCHES}', 'w') as writer:
          json.dump(batches, writer, indent=4)


  # --------------------------------------------------------------------------

  main()
  EOF
  >>>

  output {
    Array[Array[String]+]+ batches = read_json(BATCHES)
  }

  runtime {
    docker: "~{docker_image}"
    disks : "local-disk ~{disk_size} HDD"
    memory: "~{mem_size} GB"
    preemptible: preemptible
  }
}

task CheckInputWeightFiles {
    input {
        File score_weights
        Array[File] variant_weights
        String docker_image
        Int disk_size = ceil(size(score_weights, "GB") + size(variant_weights, "GB")) + 10
        Int mem_size = 2
        Int preemptible = 1
    }

    Int n_variant_weights = length(variant_weights)

    command <<<
        set -euxo pipefail

        # Make a list of PGS IDs from the variants weights files
        for file in '~{sep="' '" variant_weights}'; do
            printf -- "${file}\n" >> pgs_ids.txt
        done

        # Check for equivalent number of PGS IDs
        if [ "~{n_variant_weights}" != $(tail -n +2 "~{score_weights}" | wc -l) ]; then
            echo "ERROR: Number of PGS IDs does not match" 1>&2
            exit 1
        fi

        # Check for a matching PGS IDs in the variants weights files vs score weights file
        while read line; do
            pgs_id=$(echo "${line}" | cut -f 1)
            if [ $(grep -c ${pgs_id} pgs_ids.txt) -lt 1 ]; then
                echo "ERROR: ${pgs_id} missing from variants weights array" 1>&2
                exit 1
            fi
        done < <(tail -n +2 "~{score_weights}")
    >>>

    runtime {
        docker: docker_image
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        String input_result = "successful"
    }
}

task TrimPcaVariants {
    input {
      File    pca_variants
      File    reference
      File    query
      Boolean debug       = false
      String  docker_image
      Int     preemptible = 5
    }

    String        WORKDIR                 = "WORK"
    String        ARCHIVE                 = WORKDIR + ".tgz"
    String        OUTPUTDIR               = "OUTPUT"
    String        KEPT_PCA_VARIANTS       = OUTPUTDIR + "/kept_pca_variants.txt"
    String        WARNINGS                = OUTPUTDIR + "/WARNINGS"
    Array[String] WORKFILES               = [
                                                WORKDIR + "/PCA"
                                              , WORKDIR + "/QUERY"
                                              , WORKDIR + "/REFERENCE"
                                              , WORKDIR + "/PQ"
                                              , WORKDIR + "/PR"
                                              , WORKDIR + "/NIXQ"
                                              , WORKDIR + "/NIXR"
                                              , WORKDIR + "/NIX"
                                              , WORKDIR + "/TEMP_PCA"
                                              , WORKDIR + "/WANTED_PCA"
                                            ]
    Int           multiplier          = if debug then 20 else 10
    Int           storage             = 30 + multiplier * ceil(  size(pca_variants, "GB")
                                                               + size(query       , "GB")
                                                               + size(reference   , "GB"))

    command <<<
    set -o pipefail
    set -o errexit
    set -o nounset
    # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
    # set -o xtrace

    error() {
        local message="${1}"
        printf -- '%s\n' "${message}" >&2
        exit 1
    }

    printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
    printf -- 'INITIAL STORAGE UTILIZATION:\n'
    df --human
    printf -- '\n'

    WORKDIR='~{WORKDIR}'
    OUTPUTDIR='~{OUTPUTDIR}'
    UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
    UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf
    UNSORTEDQUERY='~{query}'         # CHROM:POS:REF:ALT from query vcf

    mkdir --verbose --parents "${WORKDIR}"
    mkdir --verbose --parents "${OUTPUTDIR}"

    WIDTH=$(( ${#WORKDIR} + 10 ))
    TEMPLATE="Creating %-${WIDTH}s ... "
    unset WIDTH

    PCA="${WORKDIR}/PCA"
    REFERENCE="${WORKDIR}/REFERENCE"
    QUERY="${WORKDIR}/QUERY"

    printf -- "${TEMPLATE}" "${PCA}"
    sort --unique "${UNSORTEDPCA}" > "${PCA}"
    printf -- 'done\n'

    printf -- "${TEMPLATE}" "${QUERY}"
    sort --unique "${UNSORTEDQUERY}" > "${QUERY}"
    printf -- 'done\n'

    printf -- "${TEMPLATE}" "${REFERENCE}"
    sort --unique "${UNSORTEDREFERENCE}" > "${REFERENCE}"
    printf -- 'done\n'

    # ---------------------------------------------------------------------------

    PQ="${WORKDIR}/PQ"
    printf -- "${TEMPLATE}" "${PQ}"
    comm -1 -2 "${QUERY}" "${PCA}" > "${PQ}"
    printf -- 'done\n'

    if [[ ! -s ${PQ} ]]
    then
        error 'No variant in the PCA variants file is mentioned in the query VCF.'
    fi

    PR="${WORKDIR}/PR"
    printf -- "${TEMPLATE}" "${PR}"
    comm -1 -2 "${REFERENCE}" "${PCA}" > "${PR}"
    printf -- 'done\n'

    if [[ ! -s ${PR} ]]
    then
        error 'No variant in the PCA variants file is mentioned in the reference VCF.'
    fi

    NPR=$( wc --lines < "${PR}" )
    NPCS=20
    if (( NPR < 2 * NPCS + 1 ))
    then

        # NB: The flashpca program, invoked by PCATasks.PerformPCA, will fail if
        #
        #     min(NVARIANTS, NSAMPLES) < 2 * NPCS + 1
        #
        # ...where NVARIANTS and NSAMPLES are the numbers of variants (rows) and
        # samples (data columns) in the reference VCF, and NPCS is the value of
        # the command's --ndim flag (currently hard-coded as 20; see
        # PCATasks.PerformPCA).  (The inequality above is derived from line 623 of
        # https://github.com/gabraham/flashpca/blob/b8044f13607a072125828547684fde8b081d6191/flashpca.cpp .)
        # In particular, flashpca will fail when NPR = NVARIANTS is less than
        # 2 * NPCS + 1.  Hence the error condition detected here.

        error 'The reference VCF mentions too few PCA variants for the PCA calculation.'
    fi

    NIXQ="${WORKDIR}/NIXQ"
    printf -- "${TEMPLATE}" "${NIXQ}"
    comm -2 -3 "${PQ}" "${PR}" > "${NIXQ}"
    printf -- 'done\n'

    NIXR="${WORKDIR}/NIXR"
    printf -- "${TEMPLATE}" "${NIXR}"
    comm -2 -3 "${PR}" "${PQ}" > "${NIXR}"
    printf -- 'done\n'

    NIX="${WORKDIR}/NIX"
    printf -- "${TEMPLATE}" "${NIX}"
    sort "${NIXQ}" "${NIXR}" > "${NIX}"
    printf -- 'done\n'

    if [[ ! -s ${NIX} ]]
    then
        exit 0
    fi

    (
        sortvariants() {
            local variants="${1}"
            sort                    \
                --field-separator=: \
                --key=1,1V          \
                --key=2,2n          \
                --key=3,3           \
                --key=4,4           \
                "${variants}"
        }

        SEPARATOR=''
        maybewarn() {
            local nixx="${1}"
            local label0="${2}"
            local label1
            if [[ ${label0} == query ]]
            then
                label1=reference
            else
                label1=query
            fi
            if [[ -s ${nixx} ]]
            then
                printf -- "${SEPARATOR}"
                if [[ -z "${SEPARATOR}" ]]
                then
                    SEPARATOR=$'\n'
                fi
                printf -- 'WARNING: '
                printf -- 'the following PCA variants from the '
                printf -- '%s VCF do not appear in the %s VCF, '  \
                          "${label0}" "${label1}"
                printf -- 'and will be removed from the '
                printf -- 'analysis.\n'
                sortvariants "${nixx}"
            fi
        }

        exec >>'~{WARNINGS}'
        maybewarn "${NIXQ}" query
        maybewarn "${NIXR}" reference
    )

    TEMP_PCA="${WORKDIR}/TEMP_PCA"
    printf -- "${TEMPLATE}" "${TEMP_PCA}"
    comm -2 -3 "${PCA}" "${NIX}" > "${TEMP_PCA}"
    printf -- 'done\n'

    WANTED_PCA="${WORKDIR}/WANTED_PCA"
    printf -- "${TEMPLATE}" "${WANTED_PCA}"
    join                                                \
        -t $'\t'                                        \
        "${TEMP_PCA}"                                   \
        <(
           perl -lpe '$_ = qq($_\t$.)' "${UNSORTEDPCA}" \
             | sort                                     \
                   --field-separator=$'\t'              \
                   --key=1,1
         )                                              \
      | sort                                            \
            --field-separator=$'\t'                     \
            --key=2,2n                                  \
      | cut --fields=1                                  \
      > "${WANTED_PCA}"
    printf -- 'done\n'

    mkdir --verbose --parents "$( dirname '~{KEPT_PCA_VARIANTS}' )"
    cp --verbose "${WANTED_PCA}" '~{KEPT_PCA_VARIANTS}'

    printf -- '\n\n## WORKDIR:\n'
    find '~{WORKDIR}' -type f | xargs ls -ltr

    (
        exec >&2
        printf -- '\n\n### LINE COUNTS:\n'
        find '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
    )

    if ~{if debug then "true" else "false"}
    then
        tar -cvzf '~{ARCHIVE}' '~{WORKDIR}'
    fi

    printf -- 'FINAL STORAGE UTILIZATION:\n'
    df --human

    if ~{if !debug then "true" else "false"}
    then
        rm -rf '~{WORKDIR}'
    fi
    >>>

    output {
      File?        kept_pca_variants    = KEPT_PCA_VARIANTS
      File?        warnings             = WARNINGS
      Array[File?] workfiles            = WORKFILES
      File?        workfiles_tgz        = ARCHIVE
    }

    runtime {
      docker: "~{docker_image}"
      disks: "local-disk ~{storage} HDD"
      preemptible: preemptible
    }
}

task TrimPcaToSubset {
  input {
    File    pca_variants
    File    reference
    File    union
    File    intersection
    Boolean debug        = false
    String  docker_image = "ubuntu:21.10"
    Int     preemptible  = 5
  }

  String        WORKDIR                 = "WORK"
  String        ARCHIVE                 = WORKDIR + ".tgz"
  String        OUTPUTDIR               = "OUTPUT"
  String        KEPT_PCA_VARIANTS       = OUTPUTDIR + "/kept_pca_variants.tsv"
  String        WARNINGS                = OUTPUTDIR + "/WARNINGS"
  Array[String] WORKFILES               = [
                                              WORKDIR + "/PCA"
                                            , WORKDIR + "/UNION"
                                            , WORKDIR + "/INTERSECTION"
                                            , WORKDIR + "/REFERENCE"
                                            , WORKDIR + "/PU"
                                            , WORKDIR + "/PI"
                                            , WORKDIR + "/PR"
                                            , WORKDIR + "/NIXQ"
                                            , WORKDIR + "/NIXR"
                                            , WORKDIR + "/NIX"
                                            , WORKDIR + "/TEMP_PCA"
                                            , WORKDIR + "/WANTED_PCA"
                                          ]
  Int           multiplier          = if debug then 20 else 10
  Int           storage             = 30 + multiplier * ceil(  size(pca_variants, "GB")
                                                             + size(union       , "GB")
                                                             + size(intersection, "GB")
                                                             + size(reference   , "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  alert() {
      local tag="${1}"
      local message="${2}"
      printf -- '%s: %s\n' "${tag}" "${message}"
  }

  error() {
      local message="${1}"
      alert 'ERROR' "${message}" >&2
      exit 1
  }

  warning() {
      local message="${1}"
      alert 'WARNING' "${message}" >> '~{WARNINGS}'
  }

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
  printf -- 'INITIAL STORAGE UTILIZATION:\n'
  df --human
  printf -- '\n'

  WORKDIR='~{WORKDIR}'
  OUTPUTDIR='~{OUTPUTDIR}'
  UNSORTEDPCA='~{pca_variants}'    # pca_variants file (no header row)
  UNSORTEDREFERENCE='~{reference}' # CHROM:POS:REF:ALT from reference vcf

  mkdir --verbose --parents "${WORKDIR}"
  mkdir --verbose --parents "${OUTPUTDIR}"

  WIDTH=$(( ${#WORKDIR} + 10 ))
  TEMPLATE="Creating %-${WIDTH}s ... "
  unset WIDTH

  PCA="${WORKDIR}/PCA"
  REFERENCE="${WORKDIR}/REFERENCE"
  UNION="${WORKDIR}/UNION"
  INTERSECTION="${WORKDIR}/INTERSECTION"

  printf -- "${TEMPLATE}" "${PCA}"
  sort --unique "${UNSORTEDPCA}" > "${PCA}"
  printf -- 'done\n'

  printf -- "${TEMPLATE}" "${REFERENCE}"
  sort --unique "${UNSORTEDREFERENCE}" > "${REFERENCE}"
  printf -- 'done\n'

  ln --symbolic --verbose '~{union}'        "${UNION}"
  ln --symbolic --verbose '~{intersection}' "${INTERSECTION}"

  PU="${WORKDIR}/PU"
  printf -- "${TEMPLATE}" "${PU}"
  comm -1 -2 "${UNION}" "${PCA}" > "${PU}"
  printf -- 'done\n'

  if [[ ! -s ${PU} ]]
  then
      error 'No variant in the PCA variants file is mentioned in the query VCFs.'
  fi

  PI="${WORKDIR}/PI"
  printf -- "${TEMPLATE}" "${PI}"
  comm -1 -2 "${INTERSECTION}" "${PCA}" > "${PI}"
  printf -- 'done\n'

  if [[ ! -s ${PI} ]]
  then
      warning 'No variant in the PCA variants file is mentioned in all the query VCFs.'
  fi

  PR="${WORKDIR}/PR"
  printf -- "${TEMPLATE}" "${PR}"
  comm -1 -2 "${REFERENCE}" "${PCA}" > "${PR}"
  printf -- 'done\n'

  if [[ ! -s ${PR} ]]
  then
      error 'No variant in the PCA variants file is mentioned in the reference VCF.'
  fi

  NPR=$( wc --lines < "${PR}" )
  NPCS=20
  if (( NPR < 2 * NPCS + 1 ))
  then

      # NB: The flashpca program, invoked by PCATasks.PerformPCA, will fail if
      #
      #     min(NVARIANTS, NSAMPLES) < 2 * NPCS + 1
      #
      # ...where NVARIANTS and NSAMPLES are the numbers of variants (rows) and
      # samples (data columns) in the reference VCF, and NPCS is the value of
      # the command's --ndim flag (currently hard-coded as 20; see
      # PCATasks.PerformPCA).  (The inequality above is derived from line 623 of
      # https://github.com/gabraham/flashpca/blob/b8044f13607a072125828547684fde8b081d6191/flashpca.cpp .)
      # In particular, flashpca will fail when NPR = NVARIANTS is less than
      # 2 * NPCS + 1.  Hence the error condition detected here.

      error 'The reference VCF mentions too few PCA variants for the PCA calculation.'
  fi

  NIXQ="${WORKDIR}/NIXQ"
  printf -- "${TEMPLATE}" "${NIXQ}"
  comm -2 -3 "${PU}" "${PR}" > "${NIXQ}"
  printf -- 'done\n'

  NIXR="${WORKDIR}/NIXR"
  printf -- "${TEMPLATE}" "${NIXR}"
  comm -2 -3 "${PR}" "${PI}" > "${NIXR}"
  # Equivalently, one could do this:
  # comm -2 -3 "${PR}" "${INTERSECTION}" > "${NIXR}"
  printf -- 'done\n'

  NIX="${WORKDIR}/NIX"
  printf -- "${TEMPLATE}" "${NIX}"
  sort "${NIXQ}" "${NIXR}" > "${NIX}"
  printf -- 'done\n'

  if [[ ! -s ${NIX} ]]
  then
      exit 0
  fi

  (
      sortvariants() {
          local variants="${1}"
          sort                    \
              --field-separator=: \
              --key=1,1V          \
              --key=2,2n          \
              --key=3,3           \
              --key=4,4           \
              "${variants}"
      }

      SEPARATOR=''
      maybewarn() {
          local nixx="${1}"
          local label0="${2}"
          local label1
          if [[ ${label0} == query ]]
          then
              label1=reference
          else
              label1=query
          fi
          if [[ -s ${nixx} ]]
          then
              printf -- "${SEPARATOR}"
              if [[ -z "${SEPARATOR}" ]]
              then
                  SEPARATOR=$'\n'
              fi
              printf -- 'WARNING: '
              printf -- 'the following PCA variants from the '
              printf -- '%s VCF do not appear in the %s VCF, '  \
                        "${label0}" "${label1}"
              printf -- 'and will be removed from the '
              printf -- 'analysis.\n'
              sortvariants "${nixx}"
          fi
      }

      exec >>'~{WARNINGS}'
      maybewarn "${NIXQ}" query
      maybewarn "${NIXR}" reference
  )

  TEMP_PCA="${WORKDIR}/TEMP_PCA"
  printf -- "${TEMPLATE}" "${TEMP_PCA}"
  comm -2 -3 "${PCA}" "${NIX}" > "${TEMP_PCA}"
  printf -- 'done\n'

  WANTED_PCA="${WORKDIR}/WANTED_PCA"
  printf -- "${TEMPLATE}" "${WANTED_PCA}"
  join                                                \
      -t $'\t'                                        \
      "${TEMP_PCA}"                                   \
      <(
         perl -lpe '$_ = qq($_\t$.)' "${UNSORTEDPCA}" \
           | sort                                     \
                 --field-separator=$'\t'              \
                 --key=1,1
       )                                              \
    | sort                                            \
          --field-separator=$'\t'                     \
          --key=2,2n                                  \
    | cut --fields=1                                  \
    > "${WANTED_PCA}"
  printf -- 'done\n'

  mkdir --verbose --parents "$( dirname '~{KEPT_PCA_VARIANTS}' )"
  cp --verbose "${WANTED_PCA}" '~{KEPT_PCA_VARIANTS}'

  printf -- '\n\n## WORKDIR:\n'
  find -L '~{WORKDIR}' -type f | xargs ls -ltr

  (
      exec >&2
      printf -- '\n\n### LINE COUNTS:\n'
      find -L '~{WORKDIR}' -type f | xargs ls -1tr | xargs wc --lines
  )

  if ~{if debug then "true" else "false"}
  then
      tar -cvzf '~{ARCHIVE}' '~{WORKDIR}'
  fi

  printf -- 'FINAL STORAGE UTILIZATION:\n'
  df --human

  if ~{if !debug then "true" else "false"}
  then
      rm -rf '~{WORKDIR}'
  fi

  >>>

  output {
    File?        kept_pca_variants    = KEPT_PCA_VARIANTS
    File?        warnings             = WARNINGS
    Array[File?] workfiles            = WORKFILES
    File?        workfiles_tgz        = ARCHIVE
    Boolean      sequencing           = true
  }

  runtime {
    docker : "~{docker_image}"
    disks : "local-disk ~{storage} HDD"
    preemptible: preemptible
  }
}
