version 1.0

workflow CramToBam {
  input {
    File    input_cram
    String  sample_name     = sub(basename(input_cram), "\\.[Cc][Rr][Aa][Mm]$", "")
    File    ref_fasta
    File    ref_fasta_index
    File    ref_dict
    String  docker          = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    Int     preemptible     = 2
  }

  call ConvertCramToBam {
    input:
        input_cram      = input_cram
      , sample_name     = sample_name
      , ref_fasta       = ref_fasta
      , ref_fasta_index = ref_fasta_index
      , ref_dict        = ref_dict
      , docker          = docker
      , preemptible     = preemptible
  }

  output {
    File output_bam = ConvertCramToBam.output_bam
    File output_bai = ConvertCramToBam.output_bai
  }
}

task ConvertCramToBam {
  input {
    File    input_cram
    String  sample_name     = sub(basename(input_cram), "\\.[Cc][Rr][Aa][Mm]$", "")
    File    ref_fasta
    File    ref_fasta_index
    File    ref_dict
    String  docker          = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    Int     preemptible     = 2
  }

  Float output_bam_size = size(input_cram, "GB") / 0.40
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int   disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20

  String output_bam_ = "~{sample_name}.bam"
  String output_bai_ = "~{sample_name}.bai"

  command <<<
  set -o errexit
  set -o pipefail

  samtools view -h -T '~{ref_fasta}' '~{input_cram}' |
  samtools view -b -o '~{output_bam_}' -
  samtools index -b '~{output_bam_}'
  mv '~{output_bam_}.bai' '~{output_bai_}'
  >>>

  runtime {
    docker      : docker
    memory      : "15 GB"
    disks       : "local-disk ~{disk_size} HDD"
    preemptible : preemptible
  }

  output {
    File output_bam = output_bam_
    File output_bai = output_bai_
  }
}
