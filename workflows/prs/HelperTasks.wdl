version 1.0

task GetBaseMemory {

  # NB: This task computes the memory (in gigabytes) required by vcf,
  # according to the recommendations given in
  # https://www.cog-genomics.org/plink/2.0/other#memory

  input {
    File vcf
  }

  Int    storage   = 20 + 2 * ceil(size(vcf, "GB"))
  String OUTPUTDIR = "OUTPUT"
  String NVARIANTS = OUTPUTDIR + "/nvariants.txt"
  String GIGABYTES = OUTPUTDIR + "/gigabytes.txt"

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  NVARIANTS="$( zgrep --count --invert-match '^#' '~{vcf}' )"

  MINGIGABYTES=8
  if (( NVARIANTS > 50000000 ))
  then
      EXCESS=$(( NVARIANTS - 50000000 ))
      DIVISOR=10000000

      # NB: The RHS in the next assignment below is equivalent to
      # ceil(EXCESS/DIVISOR); see https://stackoverflow.com/a/12536521
      EXTRA=$(( (EXCESS + DIVISOR - 1)/DIVISOR ))

      GIGABYTES=$(( MINGIGABYTES + EXTRA ))
  else
      GIGABYTES=${MINGIGABYTES}
  fi

  mkdir --verbose --parents '~{OUTPUTDIR}'
  printf -- '%d\n' "${NVARIANTS}" > '~{NVARIANTS}'
  printf -- '%d\n' "${GIGABYTES}" > '~{GIGABYTES}'
  >>>

  output {
    Int nvariants = read_int(NVARIANTS)
    Int gigabytes = read_int(GIGABYTES)
  }

  runtime {
    docker: "ubuntu:21.10"
    disks : "local-disk ~{storage} HDD"
  }
}


task RenameChromosomesInTsv {
  input {
    File    tsv
    Boolean skipheader
    File    lookup     = "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv"
  }

  Int    storage   = 20 + 2 * ceil(size(tsv, "GB"))
  String OUTPUTDIR = "OUTPUT"
  String RENAMED   = OUTPUTDIR + "/renamed_" + basename(tsv)

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
    File rename = "gs://fc-secure-9ea53c3d-d71a-4f59-92c3-63c75c622a88/reference/etc/rename_chromosomes.tsv"
  }

  Int    storage   = 20 + 2 * ceil(size(vcf, "GB"))
  String OUTPUTDIR = "OUTPUT"
  String RENAMED   = OUTPUTDIR + "/renamed_" + basename(vcf)

  command <<<
  set -o errexit
  # set -o pipefail
  # set -o nounset
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
  # set -o xtrace

  # ---------------------------------------------------------------------------

  mkdir --verbose --parents '~{OUTPUTDIR}'

  WORKDIR="$( mktemp --directory )"
  INPUTVCF="${WORKDIR}/input.vcf.gz"

  ln --symbolic --verbose '~{vcf}' "${INPUTVCF}"

  bcftools index    \
      --force       \
      --tbi         \
      "${INPUTVCF}"

  bcftools annotate            \
      --no-version             \
      --output='~{RENAMED}'    \
      --output-type=z          \
      --rename-chr='~{rename}' \
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
