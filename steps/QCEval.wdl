version 1.0

task QCEvalTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".qceval"
        String project_type
        String docker_image
        Int disk_size = ceil(size(input_vcf, "GB") * 12) + 10
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        mkdir workdir
        bcftools view --output-type v "~{input_vcf}" > workdir/input.vcf

        if [ "~{project_type}" == "WGS" ]
        then
            $MGBPMBIOFXPATH/biofx-qceval/bin/categorize_variants.py -i workdir/input.vcf -o workdir/output.vcf -g
        elif [ "~{project_type}" == "WGS_DRAGEN" ]
        then
            $MGBPMBIOFXPATH/biofx-qceval/bin/categorize_variants.py -i workdir/input.vcf -o workdir/output.vcf -d
        elif [ "~{project_type}" == "WES" ]
        then
            $MGBPMBIOFXPATH/biofx-qceval/bin/categorize_variants.py -i workdir/input.vcf -o workdir/output.vcf -e
        elif [ "~{project_type}" == "NONE" ]
        then
            $MGBPMBIOFXPATH/biofx-qceval/bin/categorize_variants.py -i workdir/input.vcf -o workdir/output.vcf -n
        else
            echo "ERROR: invalid project type ~{project_type}" 1>&2
            exit 1
        fi

        bcftools view --output-type z workdir/output.vcf > "~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
    }
}