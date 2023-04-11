version 1.0

# See DepthOfCoverage.md for detailed documentation

workflow DepthOfCoverageWorkflow {
    input {
        Boolean run_wgs = true
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        String output_basename = sub(basename(bam), "\\.(bam|BAM|cram|CRAM)$", "") + ".cov"
        File? roi_all_bed
        Array[RoiAndRefGeneFilePair?] roi_genes
        File? gene_names
        Int gatk_max_heap_gb = 31
        Int gatk_disk_size = 100
        String cov_docker_image
        String gatk_docker_image = "broadinstitute/gatk3:3.7-0"
    }

    # Only run with no interval file if so specified
    if (run_wgs) {
        call DepthOfCoverageWGSTask {
            input:
                output_basename = output_basename,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                max_heap_gb = gatk_max_heap_gb,
                disk_size = gatk_disk_size,
                docker_image = gatk_docker_image
        }
    }

    # Only run if an ROI file is provided
    if (defined(roi_all_bed)) {
        call DepthOfCoverageROITask {
            input:
                output_basename = output_basename,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                bed = select_first([roi_all_bed]),
                max_heap_gb = gatk_max_heap_gb,
                disk_size = gatk_disk_size,
                docker_image = gatk_docker_image
        }
    }

    # Only run if ROI+gene file pairs are specified
    scatter (roi_gene in select_all(roi_genes)) {
        call DepthOfCoverageGeneTask {
            input:
                output_basename = output_basename,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = bam,
                bai = bai,
                roi_bed = roi_gene.roi_bed,
                ref_gene = roi_gene.ref_gene,
                ref_gene_idx = roi_gene.ref_gene_idx,
                max_heap_gb = gatk_max_heap_gb,
                disk_size = gatk_disk_size,
                docker_image = gatk_docker_image
        }
    }

    if (defined(gene_names) && length(DepthOfCoverageGeneTask.sample_gene_summary) > 0 && defined(DepthOfCoverageROITask.sample_interval_summary)) {
        call DepthOfCoverageSummaryTask {
            input:
                sample_gene_summaries = DepthOfCoverageGeneTask.sample_gene_summary,
                sample_interval_summary = select_first([DepthOfCoverageROITask.sample_interval_summary]),
                gene_summary_file = output_basename + ".merge.sample_gene_summary.txt",
                gene_summary_unknown_file = output_basename + ".merge.sample_gene_summary.unknown.txt",
                gene_summary_entrez_file = output_basename + ".merge.sample_gene_summary.entrez.txt",
                mt_summary_file = output_basename + ".sample_mt_summary.txt",
                gene_names = select_first([gene_names]),
                docker_image = cov_docker_image
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
    File? ref_gene_idx
}

task DepthOfCoverageWGSTask {
    input {
        String output_basename
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        Int max_heap_gb
        Int disk_size
        String docker_image
    }

    command <<<
        set -euxo pipefail
        mkdir cov_out
        java -Xmx~{max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
            -I "~{bam}" \
            -ct 8 -ct 15 -ct 30 \
            -R "~{ref_fasta}" \
            -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
            -o "cov_out/~{output_basename}.nobed" \
            -omitIntervals
        ls -l cov_out
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: (max_heap_gb + 4) + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File sample_statistics = "cov_out/~{output_basename}.nobed.sample_statistics"
        File sample_summary = "cov_out/~{output_basename}.nobed.sample_summary"
    }
}

task DepthOfCoverageROITask {
    input {
        String output_basename
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        File bed
        Int max_heap_gb
        Int disk_size
        String docker_image
    }

    command <<<
        set -euxo pipefail
        mkdir cov_out
        java -Xmx~{max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
            -I "~{bam}" \
            -ct 8 -ct 15 \
            -R "~{ref_fasta}" \
            -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
            --printBaseCounts \
            -o "cov_out/~{output_basename}.roibed" \
            -L "~{bed}"
        ls -l cov_out
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: (max_heap_gb + 4) + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File sample_interval_summary = "cov_out/~{output_basename}.roibed.sample_interval_summary"
        File sample_interval_statistics = "cov_out/~{output_basename}.roibed.sample_interval_statistics"
        File sample_statistics = "cov_out/~{output_basename}.roibed.sample_statistics"
        File sample_summary = "cov_out/~{output_basename}.roibed.sample_summary"
        File sample_cumulative_coverage_counts = "cov_out/~{output_basename}.roibed.sample_cumulative_coverage_counts"
        File sample_cumulative_coverage_proportions = "cov_out/~{output_basename}.roibed.sample_cumulative_coverage_proportions"
    }
}


task DepthOfCoverageGeneTask {
    input {
        String output_basename
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        File roi_bed
        File ref_gene
        File? ref_gene_idx
        Int max_heap_gb
        Int disk_size
        String docker_image
    }

    String refg_idx = sub(basename(roi_bed), "[^0-9]", "")

    command <<<
        set -euxo pipefail

        # ensure the refseq index file name matches the expected value
        #  so that GATK does not take the time to write a new one
        if [ ! -z "~{ref_gene_idx}" -a "~{ref_gene_idx}" != "~{ref_gene}.idx" ]
        then
            cp "~{ref_gene_idx}" "~{ref_gene}.idx"
        fi
        
        mkdir cov_out
        java -Xmx~{max_heap_gb}g -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage \
            -I "~{bam}" \
            -ct 8 -ct 15 -ct 30 \
            -R "~{ref_fasta}" \
            -dt BY_SAMPLE -dcov 1000 -l INFO --omitDepthOutputAtEachBase --minBaseQuality 10 --minMappingQuality 17 --countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
            --printBaseCounts \
            -o "cov_out/~{output_basename}.refg~{refg_idx}" \
            -L "~{roi_bed}" \
            --calculateCoverageOverGenes:REFSEQ "~{ref_gene}"
        ls -l cov_out
    >>>

    runtime {
        docker: "~{docker_image}"
        memory: (max_heap_gb + 4) + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File sample_gene_summary = "cov_out/~{output_basename}.refg~{refg_idx}.sample_gene_summary"
        File sample_interval_summary = "cov_out/~{output_basename}.refg~{refg_idx}.sample_interval_summary"
        File sample_interval_statistics = "cov_out/~{output_basename}.refg~{refg_idx}.sample_interval_statistics"
        File sample_statistics = "cov_out/~{output_basename}.refg~{refg_idx}.sample_statistics"
        File sample_summary = "cov_out/~{output_basename}.refg~{refg_idx}.sample_summary"
        File sample_cumulative_coverage_counts = "cov_out/~{output_basename}.refg~{refg_idx}.sample_cumulative_coverage_counts"
        File sample_cumulative_coverage_proportions = "cov_out/~{output_basename}.refg~{refg_idx}.sample_cumulative_coverage_proportions"
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
        String docker_image
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
        docker: "~{docker_image}"
        memory: "4GB"
    }

    output {
        File gene_summary = gene_summary_file
        File gene_summary_unknown = gene_summary_unknown_file
        File gene_summary_entrez = gene_summary_entrez_file
        File mt_summary = mt_summary_file
    }
}