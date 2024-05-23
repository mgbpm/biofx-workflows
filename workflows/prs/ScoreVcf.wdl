version 1.0

workflow ScoreVcfWorkflow {
  input {
    File vcf
    File weights
  }

  call ChromosomeEncoding {
    input:
      weights = weights
  }

  call ScoreVcf {
    input:
      vcf      = vcf,
      weights  = weights,
      encoding = ChromosomeEncoding.value
  }

  output {
    File scores   = ScoreVcf.scores
    File variants = ScoreVcf.variants
  }
}

task ChromosomeEncoding {
  input {
    File weights
  }

  String outputdir = "work"
  String encoding  = outputdir + "/encoding"

  command <<<
    set -o errexit
    # set -o pipefail
    # set -o nounset
    # set -o xtrace
    # export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

    mkdir --verbose --parents '~{outputdir}'

    python3 << "EOF"
    with open('~{weights}') as reader:
        next(reader)
        chromosomes = {line.split('\t')[0].split(':')[0]
                       for line in reader}

    prefix = 'chr' if any('chr' in chromosome
                          for chromosome in chromosomes) else ''

    M  = f'{prefix}M'
    MT = f'{prefix}MT'

    code = M if M in chromosomes else MT

    with open('~{encoding}', 'w') as writer:
        writer.write(f'{code}\n')
    EOF
  >>>

  runtime {
    docker : "python:3.9.10"
  }

  output {
    String value = read_string(encoding)
  }
}


task ScoreVcf {
  input {
    File   vcf
    File   weights
    String encoding
    Int    memory   = 8
  }

  String outputdir = "work"
  String prefix    = outputdir + "/data"

  Int plink_memory = ceil(memory * 0.75 * 1000)
  Int storage      = 3 * ceil(size(vcf, "GB")) + 20

  String columns = "maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums"

  command <<<
  set -o errexit

  mkdir --verbose --parents "$( dirname '~{prefix}' )"

  /plink2                     \
      --allow-extra-chr       \
      --new-id-max-allele-len \
          1000                \
          missing             \
      --memory                \
          ~{plink_memory}     \
      --out                   \
          '~{prefix}'         \
      --output-chr            \
          '~{encoding}'       \
      --score                 \
          '~{weights}'        \
          header              \
          no-mean-imputation  \
          ignore-dup-ids      \
          list-variants       \
          cols='~{columns}'   \
      --set-all-var-ids       \
          '@:#:$1:$2'         \
      --vcf                   \
          '~{vcf}'            \
          dosage=DS
  >>>

  output {
    File scores   = "~{prefix}.sscore"
    File variants = "~{prefix}.sscore.vars"
    File log    =   "~{prefix}.log"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks:  "local-disk " + storage + " HDD"
    memory: (memory + 2) + " GB"
  }
}
