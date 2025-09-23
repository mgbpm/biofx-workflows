version 1.0

task QCEvalTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".qceval"
        String project_type
        String docker_image
        Int disk_size = ceil(size(input_vcf, "GB") * 12 + 10 +
                             (if defined(reference_fasta) then
                              size(select_first([reference_fasta]), "GB")
                              else 0) +
                             (if defined(reference_fasta_fai) then
                              size(select_first([reference_fasta_fai]), "GB")
                              else 0) +
                             (if defined(thresholds_tsv) then
                              size(select_first([thresholds_tsv]), "GB")
                              else 0) +
                             (if defined(regions_tgz) then
                              2 * size(select_first([regions_tgz]), "GB")
                              else 0))

        Int preemptible = 1
        # -------------------------------------------------------------
        # BGE_DRAGEN_TP_BINNING-specific inputs
        File? reference_fasta
        File? reference_fasta_fai
        File? thresholds_tsv
        File? regions_tgz
    }

    command <<<
        set -euxo pipefail

        mkdir workdir

        OUTPUT_VCF='~{output_basename}.vcf.gz'

        if [ "~{project_type}" == "BGE_DRAGEN_TP_BINNING" ]
        then
            REFERENCE_FASTA='workdir/reference.fasta'

            ln --symbolic --verbose '~{reference_fasta}'     "${REFERENCE_FASTA}"
            ln --symbolic --verbose '~{reference_fasta_fai}' "${REFERENCE_FASTA}.fai"

            annotate_dragen.sh       \
                "${REFERENCE_FASTA}" \
                '~{thresholds_tsv}'  \
                '~{regions_tgz}'     \
                '~{input_vcf}'       \
                "${OUTPUT_VCF}"
        else

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

            bcftools view --output-type z workdir/output.vcf > "${OUTPUT_VCF}"
        fi
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
