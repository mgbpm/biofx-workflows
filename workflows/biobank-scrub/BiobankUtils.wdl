version 1.0

# -----------------------------------------------------------------------------

task ShowEnvironment {
  input {
    String   docker_image
  }

  String    OUTPUTDIR = "OUTPUT"
  String    STDOUT    = OUTPUTDIR + "/STDOUT"

  command <<<
  set -o errexit
  set -o pipefail
  set -o nounset
  set -o xtrace
  # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

  show_environment() {
      typeset -p
      rclone version
  }

  mkdir --parents '~{OUTPUTDIR}'

  exec >&2
  show_environment | tee '~{STDOUT}'
  >>>

  output {
    File environment = STDOUT
  }

  runtime {
    docker: docker_image
  }
}

# -----------------------------------------------------------------------------
