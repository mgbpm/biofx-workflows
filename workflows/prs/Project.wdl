version 1.0

workflow ProjectWorkflow {
  input {
    File   bed
    File   bim
    File   fam
    File   loadings
    File   meansd
  }

  call Project {
    input:
      bed      = bed,
      bim      = bim,
      fam      = fam,
      loadings = loadings,
      meansd   = meansd
  }

  output {
    File projections = Project.projections
  }
}

task Project {
  input {
    File   bed
    File   bim
    File   fam
    File   loadings
    File   meansd

    Int    nthreads = 16

    Int    memory   = 8
    Int    storage  = 400
    String docker   = "us.gcr.io/broad-dsde-methods/flashpca_docker@sha256:2f3ff1614b00f9c8f271be85fd8875fbddccb7566712b537488d14a2526ccf7f"
  }

  String results = "projections.txt"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace

  # ------------------------------------------------------------------------

  check_variants() {

      local loadingsvariants="$( mktemp )"
      local meansdvariants="$( mktemp )"
      local bedvariants="$( mktemp )"

      cut --fields=1 '~{loadings}' | tail --lines=+2 > "${loadingsvariants}"
      cut --fields=1 '~{meansd}'   | tail --lines=+2 > "${meansdvariants}"
      cut --fields=2 '~{bim}'                        > "${bedvariants}"

      if ! diff --brief "${loadingsvariants}" "${meansdvariants}" >/dev/null
      then
          cat <<EOERROR >&2
ERROR: Loadings and meansd files do not contain the same variants (or they are
       not in the same order)
EOERROR
          return 1
      fi

      if ! diff --brief "${meansdvariants}" "${bedvariants}" >/dev/null
      then
          cat <<EOERROR >&2
ERROR: The variants in the .bim differ from those in the PCA files (or they are
       not in the same order)
EOERROR
          return 1
      fi
  }

  # ------------------------------------------------------------------------

  check_variants

  WORKDIR="${PWD}/work"
  BEDPREFIX="${WORKDIR}/data"

  mkdir --verbose --parents "${WORKDIR}"

  # The following cp commands ensure that all the bed-related files have the
  # same prefix and the correct extensions.
  cp '~{bed}' "${BEDPREFIX}.bed"
  cp '~{bim}' "${BEDPREFIX}.bim"
  cp '~{fam}' "${BEDPREFIX}.fam"

  # ------------------------------------------------------------------------

  ~/flashpca/flashpca           \
    --bfile      "${BEDPREFIX}" \
    --inload     '~{loadings}'  \
    --inmeansd   '~{meansd}'    \
    --numthreads ~{nthreads}    \
    --outproj    '~{results}'   \
    --project                   \
    --verbose
  >>>

  output {
    File projections = results
  }

  runtime {
    docker: docker
    disks: "local-disk " + storage + " HDD"
    memory: memory + " GB"
  }
}
