version 1.0

# See DepthOfCoverage.md for detailed documentation

workflow DepthOfCoverageWorkflow {
    input {
        Boolean run_wgs = true
        String sample_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        File? roi_all_bed
        Array[RoiAndRefGeneFilePair?] roi_genes
        File? gene_names
        Int gatk_max_heap_gb = 31
        Int gatk_disk_size_gb = 100
        Int gatk_num_cpu = floor(gatk_max_heap_gb / 4)
        String mgbpmbiofx_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/base:latest"
    }

    # Only run with no interval file if so specified
    if (run_wgs) {
        call DepthOfCoverageWGSTask {
            input:
                sample_name = sample_name,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                gatk_max_heap_gb = gatk_max_heap_gb,
                gatk_disk_size_gb = gatk_disk_size_gb,
                gatk_num_cpu = gatk_num_cpu
        }
    }

    # Only run if an ROI file is provided
    if (defined(roi_all_bed)) {
        call DepthOfCoverageROITask {
            input:
                sample_name = sample_name,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                bed = select_first([roi_all_bed]),
                gatk_max_heap_gb = gatk_max_heap_gb,
                gatk_disk_size_gb = gatk_disk_size_gb,
                gatk_num_cpu = gatk_num_cpu
        }
    }

    # Only run if ROI+gene file pairs are specified
    scatter (roi_gene in select_all(roi_genes)) {
        call DepthOfCoverageGeneTask {
            input:
                sample_name = sample_name,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                roi_bed = roi_gene.roi_bed,
                ref_gene = roi_gene.ref_gene,
                gatk_max_heap_gb = gatk_max_heap_gb,
                gatk_disk_size_gb = gatk_disk_size_gb,
                gatk_num_cpu = gatk_num_cpu
        }
    }

    if (defined(gene_names) && length(DepthOfCoverageGeneTask.sample_gene_summary) > 0 && defined(DepthOfCoverageROITask.sample_interval_summary)) {
        call DepthOfCoverageSummaryTask {
            input:
                sample_gene_summaries = DepthOfCoverageGeneTask.sample_gene_summary,
                sample_interval_summary = select_first([DepthOfCoverageROITask.sample_interval_summary]),
                gene_summary_file = sample_name + ".cov.merge.sample_gene_summary.txt",
                gene_summary_unknown_file = sample_name + ".cov.merge.sample_gene_summary.unknown.txt",
                gene_summary_entrez_file = sample_name + ".cov.merge.sample_gene_summary.entrez.txt",
                mt_summary_file = sample_name + ".cov.sample_mt_summary.txt",
                gene_names = select_first([gene_names]),
                mgbpmbiofx_docker_image = mgbpmbiofx_docker_image
        }
    }

    output {
        File? wgs_sample_summary = DepthOfCoverageWGSTask.sample_summary
        File? wgs_sample_statistics = DepthOfCoverageWGSTask.sample_statistics
        File? roi_sample_interval_summary = DepthOfCoverageROITask.sample_interval_summary
        File? roi_sample_interval_statistics = DepthOfCoverageROITask.sample_interval_statistics
        File? roi_sample_statistics = DepthOfCoverageROITask.sample_statistics
        File? roi_sample_summary = DepthOfCoverageROITask.sample_summary
        File? roi_sample_cumulative_coverage_counts = DepthOfCoverageROITask.sample_cumulative_coverage_counts
        File? roi_sample_cumulative_coverage_proportions = DepthOfCoverageROITask.sample_cumulative_coverage_proportions
        File? mt_summary = DepthOfCoverageSummaryTask.mt_summary
        File? gene_summary = DepthOfCoverageSummaryTask.gene_summary
        File? gene_summary_unknown = DepthOfCoverageSummaryTask.gene_summary_unknown
        File? gene_summary_entrez = DepthOfCoverageSummaryTask.gene_summary_entrez
    }
}

struct RoiAndRefGeneFilePair {
    File roi_bed
    File ref_gene
}

task DepthOfCoverageWGSTask {
    input {
        String sample_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        Int gatk_max_heap_gb
        Int gatk_disk_size_gb
        Int gatk_num_cpu
    }

    command <<<
        set -euxo pipefail
        mkdir cov_out
        java -Xmx~{gatk_max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
            -I "~{bam}" \
            -ct 8 -ct 15 -ct 30 \
            -R "~{ref_fasta}" \
            -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
            -o "cov_out/~{sample_name}.cov.nobed" \
            -omitIntervals \
            -nt ~{gatk_num_cpu}
        ls -l cov_out
    >>>

    runtime {
        docker: "broadinstitute/gatk3:3.7-0"
        memory: "~{gatk_max_heap_gb + 4}GB"
        cpu: gatk_num_cpu
        disks: "local-disk ~{gatk_disk_size_gb} SSD"
    }

    output {
        File sample_statistics = "cov_out/~{sample_name}.cov.nobed.sample_statistics"
        File sample_summary = "cov_out/~{sample_name}.cov.nobed.sample_summary"
    }
}

task DepthOfCoverageROITask {
    input {
        String sample_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        File bed
        Int gatk_max_heap_gb
        Int gatk_disk_size_gb
        Int gatk_num_cpu
    }

    command <<<
        set -euxo pipefail
        mkdir cov_out
        java -Xmx~{gatk_max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
            -I "~{bam}" \
            -ct 8 -ct 15 \
            -R "~{ref_fasta}" \
            -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
            --printBaseCounts \
            -o "cov_out/~{sample_name}.cov.roibed" \
            -L "~{bed}" \
            -nt ~{gatk_num_cpu}
        ls -l cov_out
    >>>

    runtime {
        docker: "broadinstitute/gatk3:3.7-0"
        memory: "~{gatk_max_heap_gb + 4}GB"
        cpu: gatk_num_cpu
        disks: "local-disk ~{gatk_disk_size_gb} SSD"
    }

    output {
        File sample_interval_summary = "cov_out/~{sample_name}.cov.roibed.sample_interval_summary"
        File sample_interval_statistics = "cov_out/~{sample_name}.cov.roibed.sample_interval_statistics"
        File sample_statistics = "cov_out/~{sample_name}.cov.roibed.sample_statistics"
        File sample_summary = "cov_out/~{sample_name}.cov.roibed.sample_summary"
        File sample_cumulative_coverage_counts = "cov_out/~{sample_name}.cov.roibed.sample_cumulative_coverage_counts"
        File sample_cumulative_coverage_proportions = "cov_out/~{sample_name}.cov.roibed.sample_cumulative_coverage_proportions"
    }
}


task DepthOfCoverageGeneTask {
    input {
        String sample_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        File roi_bed
        File ref_gene
        Int gatk_max_heap_gb
        Int gatk_disk_size_gb
        Int gatk_num_cpu
    }

    String refg_idx = sub(basename(roi_bed), "[^0-9]", "")

    command <<<
        set -euxo pipefail
        mkdir cov_out
        java -Xmx~{gatk_max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
            -I "~{bam}" \
            -ct 8 -ct 15 -ct 30 \
            -R "~{ref_fasta}" \
            -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
            --printBaseCounts \
            -o "cov_out/~{sample_name}.cov.refg~{refg_idx}" \
            -L "~{roi_bed}" \
            --calculateCoverageOverGenes:REFSEQ "~{ref_gene}" \
            -nt ~{gatk_num_cpu}
        ls -l cov_out
    >>>

    runtime {
        docker: "broadinstitute/gatk3:3.7-0"
        memory: "~{gatk_max_heap_gb + 4}GB"
        cpu: gatk_num_cpu
        disks: "local-disk ~{gatk_disk_size_gb} SSD"
    }

    output {
        File sample_gene_summary = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_gene_summary"
        File sample_interval_summary = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_interval_summary"
        File sample_interval_statistics = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_interval_statistics"
        File sample_statistics = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_statistics"
        File sample_summary = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_summary"
        File sample_cumulative_coverage_counts = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_cumulative_coverage_counts"
        File sample_cumulative_coverage_proportions = "cov_out/~{sample_name}.cov.refg~{refg_idx}.sample_cumulative_coverage_proportions"
    }
}

task DepthOfCoverageSummaryTask {
    input {
        Array[File] sample_gene_summaries
        File sample_interval_summary
        File gene_names
        String gene_summary_file
        String gene_summary_unknown_file
        String gene_summary_entrez_file
        String mt_summary_file
        String mgbpmbiofx_docker_image = "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/base:latest"
    }

    command <<<
        set -euxo pipefail
        # Split the gene summaries into known and unknown genes
        head -1 "~{sample_gene_summaries[0]}" > "~{gene_summary_file}"
        for file in ~{sep=' ' sample_gene_summaries}; do 
            tail -n +2 "$file" | grep -v UNKNOWN >> "~{gene_summary_file}"
        done
        cat ~{sep=' ' sample_gene_summaries} | grep UNKNOWN > "~{gene_summary_unknown_file}"

        # Enrich the gene summary file with Entrez IDs
        $MGBPMBIOFXPATH/biofx-depth-of-coverage/bin/add_entrez.py -i "~{gene_summary_file}" -o "~{gene_summary_entrez_file}" -m "~{gene_names}"

        # Extract MT interval from summary file to a MT summary file
        egrep "^(MT|chrM)" "~{sample_interval_summary}" >>/dev/null 2>&1
        if [ $? -eq 0 ]
        then
            head -1 "~{sample_interval_summary}" > "~{mt_summary_file}"
            grep -E "^(MT|chrM)" "~{sample_interval_summary}" >> "~{mt_summary_file}"
        fi
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "4GB"
    }

    output {
        File gene_summary = gene_summary_file
        File gene_summary_unknown = gene_summary_unknown_file
        File gene_summary_entrez = gene_summary_entrez_file
        File mt_summary = mt_summary_file
    }
}