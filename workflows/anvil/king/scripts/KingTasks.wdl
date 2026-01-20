version 1.0

task FilterVcfTask {
    input {
        File input_vcf
        File? input_bed
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String docker_image
        Int addldisk = 10
        Int preemptible = 1
    }

    Int vcf_size = ceil((size(input_vcf, "GB") * 2))
    Int bed_size = ceil((size(input_bed, "GB")))
    Int final_disk_size = vcf_size + bed_size + addldisk

    command <<<
        set -euxo pipefail

        if [ -f "~{input_bed}" ]; then
            bcftools filter                    \
                --regions-file "~{input_bed}"  \
                --output-type z "~{input_vcf}" \
                --output "~{output_basename}.vcf.gz"

            num_snps=$(wc -l < "~{input_bed}")
            num_samples=$(bcftools query -l "~{output_basename}.vcf.gz" | wc -l)
        else
            num_snps=$(zgrep -v "#" "~{input_vcf}" | wc -l)
            num_samples=$(bcftools query -l "~{input_vcf}" | wc -l)
        fi

        printf -- ${num_snps} > "num_snps.txt"
        printf -- ${num_samples} > "num_samples.txt"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + final_disk_size + " SSD"
        preemptible: preemptible
    }

    output {
        File? output_vcf_gz = "~{output_basename}.vcf.gz"
        Int num_snps = read_int("num_snps.txt")
        Int num_samples = read_int("num_samples.txt")
    }

}

task MergeVcfsTask {
    input {
        Array[File] input_vcfs
        Array[File?] input_vcfs_idx
        File? input_bed
        String output_basename
        String docker_image
        Int addldisk = 10
        Int mem_size = if size(input_vcfs, "GB") > 10 then 8 else 4 
        Int preemptible = 1
    }

    Int vcf_size = ceil(size(input_vcfs, "GB") * 3)
    Int idx_size = ceil(size(input_vcfs_idx, "GB"))
    Int bed_size = ceil(size(input_bed, "GB"))
    Int final_disk_size = vcf_size + idx_size + bed_size + addldisk

    command <<<
        set -euxo pipefail

        for c in '~{sep="' '" input_vcfs}'
        do
            echo $c >> "merge_these_files.txt"
        done

        if [ -f "~{input_bed}" ]; then
            bcftools merge                          \
                --merge none                        \
                --file-list "merge_these_files.txt" \
                --regions-file "~{input_bed}"       \
                --output-type z                     \
                --output "~{output_basename}.vcf.gz"

            num_snps=$(wc -l < "~{input_bed}")
        else
            bcftools merge                          \
                --merge none                        \
                --file-list "merge_these_files.txt" \
                --output-type z                     \
                --output "~{output_basename}.vcf.gz"
            
            num_snps=$(zgrep -v "#" "~{output_basename}.vcf.gz" | wc -l)
        fi

        printf -- ${num_snps} > "num_snps.txt"
        num_samples=$(bcftools query -l "~{output_basename}.vcf.gz" | wc -l)
        printf -- ${num_samples} > "num_samples.txt"
    >>>

    runtime {
        docker: docker_image
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File merged_vcf = "~{output_basename}.vcf.gz"
        Int num_snps = read_int("num_snps.txt")
        Int num_samples = read_int("num_samples.txt")
    }
}

task Vcf2BedTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "")
        String docker_image
        Int addldisk = 10
        Int plink_mem = 4
        Int preemptible = 1
    }

    # PLINK mem should be in MB
    Int actual_plink_mem = plink_mem * 1000
    # mem_size is in GB -- by default, PLINK would be half of alloted mem_size
    Int mem_size = plink_mem * 2
    # Estimate disk space from VCFs
    Int final_disk_size = ceil(size(input_vcf, "GB") * 4) + addldisk

    command <<<
        plink --vcf "~{input_vcf}"  \
            --make-bed              \
            --memory "~{actual_plink_mem}" \
            --out "~{output_basename}"
    >>>

    runtime {
        docker: docker_image
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File bed_file = "~{output_basename}.bed"
        File bim_file = "~{output_basename}.bim"
        File fam_file = "~{output_basename}.fam"
    }
}

task RunKingTask {
    input {
        File bed_file
        File fam_file
        File bim_file
        Int degree = 3
        String flag
        String output_basename
        String docker_image
        Int addldisk = 10
        Int cpu = 2
        Int mem_size = 4
        Int preemptible = 2
    }

    Int bed_size = ceil(size(bed_file, "GB") * 2)
    Int bim_size = ceil(size(bim_file, "GB"))
    Int fam_size = ceil(size(fam_file, "GB"))
    Int final_disk_size = bed_size + bim_size + fam_size + addldisk

    command <<<
        set -euxo pipefail

        if ["~{flag}" == "duplicate"]; then
            king --"~{flag}"                  \
                -b "~{bed_file}"              \
                --bim "~{bim_file}"           \
                --fam "~{fam_file}"           \
                --cpus "~{cpu}"               \
                --prefix "~{output_basename}"
        else
            king --"~{flag}"                  \
                -b "~{bed_file}"              \
                --bim "~{bim_file}"           \
                --fam "~{fam_file}"           \
                --cpus "~{cpu}"               \
                --prefix "~{output_basename}" \
                --degree "~{degree}"
        fi

    >>>

    runtime {
        cpu: cpu
        docker: docker_image
        disks: "local-disk " + final_disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File? seg_output = "~{output_basename}.seg"
        File? con_output = "~{output_basename}.con"
        File? kin_output = "~{output_basename}.kin"
        File? kin0_output = "~{output_basename}.kin0"
    }
}
