version 1.0

task SortVCFTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".sorted"
        Boolean output_index = false
        String index_format = "csi" # either "csi" or "tbi"
        File? empty_output_placeholder
        String docker_image
        # bcftools uses uncompressed temp files for sorting
        #  so disk needs to be able to hold uncompressed VCF data
        Int disk_size = ceil(size(input_vcf, "GB") * 23) + 10
        Int mem_size = if size(input_vcf, "GB") > 10 then 8 else 4 
        Float sort_max_mem = (if size(input_vcf, "GB") > 10 then 8 else 4) * 0.75
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        bcftools sort --max-mem ~{sort_max_mem}G --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"

        # Build a new index file if desired
        if [ "~{output_index}" == "true" ]
        then
            bcftools index --~{index_format} "~{output_basename}.vcf.gz"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
        File? output_vcf_idx = if output_index then "~{output_basename}.vcf.gz.~{index_format}" else empty_output_placeholder
    }
}

task FilterVCFWithBEDTask {
    input {
        File input_vcf
        File? input_vcf_idx
        File input_bed
        String regions_or_targets = ""
        String target_overlap = "record"
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".filtered"
        Boolean output_index = false
        String index_format = "csi" # either "csi" or "tbi"
        String docker_image
        Int disk_size = ceil((size(input_vcf, "GB") * 2) + size(input_bed, "GB")) + 10
        Int preemptible = 1
        File? empty_output_placeholder
    }

    command <<<
        set -euxo pipefail

        # if VCF index was provided, make sure it is resolvable from the VCF file path
        if [[ -n "~{input_vcf_idx}" && ! -f "~{input_vcf}.csi" && ! -f "~{input_vcf}.tbi" ]]
        then
            if [[ "~{input_vcf_idx}" =~ \.tbi$ ]]
            then
                ln -s "~{input_vcf_idx}" "~{input_vcf}.tbi"
            else
                ln -s "~{input_vcf_idx}" "~{input_vcf}.csi"
            fi
        fi

        # Determine whether to use targets or regions filtering option
        filter_type="~{regions_or_targets}"
        if [[ $filter_type != "targets" && $filter_type != "regions" ]]
        then
            # If VCF index was provided, determine how many regions there
            #   are for the chromosomes in VCF file
            if [ -n "~{input_vcf_idx}" ]
            then
                region_count=$((0))
                for seq in $(bcftools index -s "~{input_vcf}" | cut -f1)
                do
                    seq_region_count=$(( $(grep -c "^$seq\\b" "~{input_bed}") ))
                    region_count=$((region_count + seq_region_count))
                done
            # Else just count the total number of regions in the bed file
            else
                region_count=$(( $(wc -l < "~{input_bed}") ))
            fi

            # Less than 4500 regions, use regions, otherwise targets
            if [ $region_count -lt 4500 ]
            then
                filter_type="regions"
            else
                filter_type="targets"
            fi
        fi

        # If filtering by regions, but no VCF index, create one
        if [[ $filter_type == "regions" && ! -f "~{input_vcf_idx}" ]]
        then
            bcftools index "~{input_vcf}"
        fi

        # Filter the vcf
        bcftools filter --${filter_type}-file "~{input_bed}" --${filter_type}-overlap ~{target_overlap} --output-type z "~{input_vcf}" > "~{output_basename}.vcf.gz"

        # Build a new index file if desired
        if [ "~{output_index}" == "true" ]
        then
            bcftools index --~{index_format} "~{output_basename}.vcf.gz"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: preemptible
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
        File? output_vcf_idx = if output_index then "~{output_basename}.vcf.gz.~{index_format}" else empty_output_placeholder
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

task ConcatVCFsTask {
    input {
        Array[File] input_vcfs
        Boolean sorted = false
        Boolean output_index = false
        String index_format = "csi" # either "csi" or "tbi"
        String output_basename
        String docker_image
        # bcftools uses uncompressed temp files for sorting
        #  so disk needs to be able to hold uncompressed VCF data
        Int disk_size = ceil(size(input_vcfs, "GB") * 23) + 10
        Int mem_size = if size(input_vcfs, "GB") > 10 then 8 else 4 
        Float sort_max_mem = (if size(input_vcfs, "GB") > 10 then 8 else 4) * 0.75
        Int preemptible = 1
        File? empty_output_placeholder
    }

    command <<<
        set -euxo pipefail

        # Show vm stats
        bash -c 'while [ 1 -eq 1 ]; do vmstat 1 2; df -H; sleep 300; done' &

        # Create list of files to concat
        for c in '~{sep="' '" input_vcfs}'
        do
            echo $c >> "concat_these_files.txt"
        done

        # Concat files in list and sort if desired
        if [ "~{sorted}" == "true" ]
        then
            bcftools concat --file-list "concat_these_files.txt" --output-type z > "~{output_basename}_tmp.vcf.gz"

            bcftools sort --max-mem ~{sort_max_mem}G --output-type z "~{output_basename}_tmp.vcf.gz" > "~{output_basename}.vcf.gz"
        else
            bcftools concat --file-list "concat_these_files.txt" --output-type z > "~{output_basename}.vcf.gz"
        fi

        # Build a new index file if desired
        if [ "~{output_index}" == "true" ]
        then
            bcftools index --~{index_format} "~{output_basename}.vcf.gz"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
        File? output_vcf_idx = if output_index then "~{output_basename}.vcf.gz.~{index_format}" else empty_output_placeholder
    }
}

task ExtractSamplesFromVCFTask {
    input {
        File input_vcf
        Array[String] sample_ids
        String min_ac = "1"
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".subset"
        String docker_image
        Int disk_size = 10 + ceil(size(input_vcf, "GB") * 2)
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # subset the vcf to the samples and create a file for the results
        bcftools view "~{input_vcf}" --min-ac ~{min_ac} --samples "~{sep=',' sample_ids}" --output-type z > "~{output_basename}.vcf.gz"
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: preemptible
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
    }
}

task ConvertBCFTask {
    input {
        File input_bcf
        Boolean output_index = false
        String index_format = "csi" # either "csi" or "tbi"
        String output_basename = sub(basename(input_bcf), "\\.(bcf|BCF|bcf.gz|BCF.GZ|bcf.bgz|BCF.BGZ)$", "")
        String docker_image
        Int disk_size = 10 + ceil(size(input_bcf, "GB") * 2.5)
        Int preemptible = 1
        File? empty_output_placeholder
    }

    command <<<
        set -euxo pipefail

        # Convert based on desired file output format and create index if desired
        bcftools convert --output-type z --output "~{output_basename}.vcf.gz" "~{input_bcf}"
        if [ "~{output_index}" == "true" ]
        then
            bcftools index --~{index_format} "~{output_basename}.vcf.gz"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        File? output_vcf_gz = "~{output_basename}.vcf.gz"
        File? output_vcf_idx = if output_index then "~{output_basename}.vcf.gz.~{index_format}" else empty_output_placeholder
    }
}

task MakeCollectiveVCFTask {
    input {
        File input_vcf
        String output_basename = sub(basename(input_vcf), "\\.(vcf|VCF|vcf.gz|VCF.GZ|vcf.bgz|VCF.BGZ)$", "") + ".collective"
        String collective_sample_name = "CollectiveSample"
        String collective_gt_call = "0/0"
        String docker_image
        Int disk_size = ceil(size(input_vcf, "GB") * 1.5) + 10
        Int preemptible = 1
        Int mem_size = 4
    }

    command <<<
        pip install vcfpy

        # For each record, replace all samples with a fake sample
        python -c 'import vcfpy
        from collections import OrderedDict
        with vcfpy.Reader.from_path("~{input_vcf}") as rd:
            rd.header.parsed_samples=set([])
            outhdr = rd.header.copy()
            outhdr.samples = vcfpy.SamplesInfos(["~{collective_sample_name}"])
            with vcfpy.Writer.from_path("~{output_basename}.vcf.gz", outhdr) as wr:
                for rec in rd:
                    rec.FORMAT = ["GT"]
                    rec.update_calls([vcfpy.Call("~{collective_sample_name}", OrderedDict([("GT", "~{collective_gt_call}")]))])
                    wr.write_record(rec)'
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_size + "GB"
        preemptible: preemptible
    }

    output {
        File output_vcf_gz = "~{output_basename}.vcf.gz"
    }
}