version 1.0
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# -----------------------------------------------------------------------------

workflow Main {
  input {
    String image
    File   script
  }

  call RunScript {
    input:
      image  = image
    , script = script
  }

  output {
    File result = RunScript.result
  }
}

task RunScript {
  input {
    String image
    File   script
  }

  String result = "work/result"

  command <<<
  chmod +x '~{script}'

  mkdir --parents "$( dirname '~{result}' )"

  '~{script}' > '~{result}'
  >>>

  output {
    File result = result
  }

  runtime {
    docker: image
  }
}
