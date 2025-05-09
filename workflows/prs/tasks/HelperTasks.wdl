version 1.0

task GetBaseMemory {

  # NB: This task computes the memory (in gigabytes) required by vcf,
  # according to the recommendations given in
  # https://www.cog-genomics.org/plink/2.0/other#memory

  input {
    File? vcf
    Int?  nvariants
  }

  Int     storage   = 20 + 2 * ceil(size(vcf, "GB"))
  Boolean ERROR     = defined(vcf) == defined(nvariants)
  String  OUTPUTDIR = "OUTPUT"
  String  NVARIANTS = OUTPUTDIR + "/nvariants.txt"
  String  GIGABYTES = OUTPUTDIR + "/gigabytes.txt"

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
    disks : "local-disk ~{storage} HDD"
    docker: "python:3.11"
  }
}


task RenameChromosomesInTsv {
  input {
    File    tsv
    Boolean skipheader
    File    lookup     = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
  }

  Int    storage = 20 + 2 * ceil(size(tsv, "GB"))
  String RENAMED = "OUTPUT/renamed/" + basename(tsv)

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
    docker: "python:3.11"
    disks : "local-disk ~{storage} HDD"
  }

}


task RenameChromosomesInVcf {
  input {
    File vcf
    File rename = "gs://lmm-reference-data/prsmix/reference/rename_chromosomes.tsv"
  }

  Int    storage = 20 + 2 * ceil(size(vcf, "GB"))
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
      --rename-chr='~{rename}'         \
      --set-id '%CHROM:%POS:%REF:%ALT' \
      "${INPUTVCF}"

  >>>

  output {
    File renamed = RENAMED
  }

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks : "local-disk ~{storage} HDD"
  }
}

task SubsetVcf {
  input {
    File    inputvcf
    File    regions
    String  label     = "data"
    Boolean nocleanup = false
  }

  String OUTPUTDIR = "OUTPUT"
  String OUTPUTVCF = OUTPUTDIR + "/" + label + ".vcf.gz"
  String NREGIONS  = OUTPUTDIR + "/NREGIONS"
  Int    storage   = 20 + 3 * ceil(size(inputvcf, "GB"))

  command <<<
  set -o pipefail
  set -o errexit
  set -o nounset
  export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  set -o xtrace

  # ---------------------------------------------------------------------------

  printf -- 'SPECIFIED STORAGE: %d GB\n\n' '~{storage}'
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
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks : "local-disk ~{storage} HDD"
  }
}

task Union {
  input {
    Array[File]+ lists
    Int          storage
  }

  String OUTPUTDIR = "OUTPUT"
  String RESULT    = OUTPUTDIR + "/RESULT"

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
    docker: "ubuntu:21.10"
    disks : "local-disk " + (20 + storage) + " HDD"
  }
}

task Intersection {
  input {
    Array[File]+ lists
    Int          storage
  }

  String OUTPUTDIR = "OUTPUT"
  String RESULT    = OUTPUTDIR + "/RESULT"

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
    docker: "ubuntu:21.10"
    disks : "local-disk " + (20 + storage) + " HDD"
  }
}

task ListShards {
  input {
    String source
    String workspace
    String docker_image
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
    docker: docker_image
  }
}

task MakeBatches {
  input {
    File          cases
    Int           nbatches
    Array[String] exclude   = []
    Boolean       noshuffle = false
    Int?          seed
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
    docker: "python:3.9.10"
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