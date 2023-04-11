version 1.0

task SortVCFTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".sorted"
        String docker_image
        Int disk_size = 20
    }

    command <<<
        bcftools sort --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + "GB HDD"
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
    }
}

task FilterVCFWithBEDTask {
    input {
        File input_vcf
        File input_bed
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".filtered"
        String docker_image
        Int disk_size = 20
    }

    command <<<
        bcftools filter --targets-file "~{input_bed}" --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + "GB HDD"
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
        Int disk_size = 50
    }

    command <<<
        if "~{annotations_file}.idx" != "~{annotations_idx_file}" ]
        then
            mv "~{annotations_idx_file}" "~{annotations_file}.tbi"
        fi
        bcftools annotate -a "~{annotations_file}" -h "~{headers_file}" -c "~{column_list}" --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + "GB HDD"
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
    }
}