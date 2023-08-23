version 1.0

task SortVCFTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".sorted"
        String docker_image
        Int disk_size = ceil(size(input_vcf, "GB") * 2) + 10
        Int preemptible = 1
    }

    command <<<
        bcftools sort --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"
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

task FilterVCFWithBEDTask {
    input {
        File input_vcf
        File input_bed
        String target_overlap = "record"
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".filtered"
        String docker_image
        Int disk_size = ceil((size(input_vcf, "GB") * 2) + size(input_bed, "GB")) + 10
        Int preemptible = 1
    }

    command <<<
        bcftools filter --targets-file "~{input_bed}" --targets-overlap ~{target_overlap} --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"
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

task AnnotateVCFTask {
    input {
        File input_vcf
        File annotations_file
        File annotations_idx_file
        File headers_file
        String column_list
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".annotated"
        String docker_image
        Int disk_size = ceil((size(input_vcf, "GB") * 2.5) + size(annotations_file, "GB") + size(annotations_idx_file, "GB")) + 10
        Int preemptible = 1
    }

    command <<<
        if [ "~{annotations_file}.tbi" != "~{annotations_idx_file}" ]
        then
            ln -s "~{annotations_idx_file}" "~{annotations_file}.tbi"
        fi
        bcftools annotate -a "~{annotations_file}" -h "~{headers_file}" -c "~{column_list}" --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"
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