version 1.0

workflow cram_bam_to_fastq_workflow {
    input {
        File input_file
        File? refFasta
        File? refFastaIndex
        Int runtime_cpu = 4
        Int runtime_memory = 16
        Int runtime_preemptible = 3
        Float runtime_disk_multiplier = 12
    }
    
    # Determine file type based on extension
    String file_extension = sub(basename(input_file), ".*\\.", "")
    Boolean is_cram = file_extension == "cram"
    
    # Convert CRAM to BAM if input is CRAM
    if (is_cram) {
        # Call project specific cram_to_bam task
        call CBTFQ_cram_to_bam {
            input:
                cram = input_file,
                refFasta = select_first([refFasta]),
                refFastaIndex = select_first([refFastaIndex]),
                runtime_cpu = runtime_cpu,
                runtime_memory = runtime_memory,
                runtime_preemptible = runtime_preemptible,
                runtime_disk_multiplier = runtime_disk_multiplier
        }
    }
    
    # Use converted BAM or original BAM file
    File bam_file = select_first([cram_to_bam.bam_output, input_file])
    
    # Convert BAM to FASTQ
    call bam_to_fastq {
        input:
            bam = bam_file,
            runtime_cpu = runtime_cpu,
            runtime_memory = runtime_memory,
            runtime_preemptible = runtime_preemptible,
            runtime_disk_multiplier = runtime_disk_multiplier
    }
    
    output {
        File flagstat = bam_to_fastq.flagstat
        File read1 = bam_to_fastq.read1
        File read2 = bam_to_fastq.read2
        File read0 = bam_to_fastq.read0
        File singleton_reads = bam_to_fastq.singleton_reads
        Int count_0 = bam_to_fastq.count_0
        Int count_single = bam_to_fastq.count_single
        String input_type = if is_cram then "CRAM" else "BAM"
    }
}

task CBTFQ_cram_to_bam {
    input {
        File cram
        File refFasta
        File refFastaIndex
        Int runtime_cpu
        Int runtime_memory
        Int runtime_preemptible
        Float runtime_disk_multiplier
    }
    
    String name = sub(basename(cram), '\\.cram$', '')
    Int runtime_disk_gb = ceil(size(cram, 'GB') * runtime_disk_multiplier)
    
    command <<<
        set -euxo pipefail
        
        # Convert CRAM to BAM
        samtools view -b -T ~{refFasta} -@ ~{runtime_cpu} -o ~{name}.bam ~{cram}
        
        # Index the BAM file
        samtools index ~{name}.bam
    >>>
    
    output {
        File bam_output = '${name}.bam'
        File bam_index = '${name}.bam.bai'
    }
    
    runtime {
        disks: 'local-disk ${runtime_disk_gb} HDD'
        cpu: '${runtime_cpu}'
        memory: '${runtime_memory} GB'
        docker: 'quay.io/biocontainers/samtools:1.17--hd87286a_2'
        preemptible: '${runtime_preemptible}'
    }
}

task bam_to_fastq {
    input {
        File bam
        Int runtime_cpu
        Int runtime_memory
        Int runtime_preemptible
        Float runtime_disk_multiplier
    }
    
    String name = sub(basename(bam), '\\.bam$', '')
    Int collate_threads = runtime_cpu
    Int runtime_disk_gb = ceil(size(bam, 'GB') * runtime_disk_multiplier)
    
    command <<<
        set -euxo pipefail
        
        # Collect metrics
        samtools flagstat -@ ~{runtime_cpu} ~{bam} > ~{name}_flagstat.txt
        
        #Samtools collate Flags
        #Purpose: Reorganizes BAM records to group read pairs together
        #-O: Output to stdout (pipes to next command)
        #-u: Output uncompressed BAM (faster for piping)
        #-n 128: Use 128 as the number of reads to buffer in memory
        #--threads: Number of threads for parallel processing
        #temp: Temporary filename prefix for intermediate files
        #~{bam}: Input BAM file

        #Samtools fastq Flags
        #Purpose: Converts the collated BAM stream into FASTQ format, handling different read types.
        #Core flags
        #-F 0x900: Filter out unwanted reads
        #    0x900 = 0x100 | 0x800 (hexadecimal)
        #    0x100 = secondary alignments
        #    0x800 = supplementary alignments
        #    This keeps only primary alignments, avoiding duplicates
        #-n: Don't append /1 and /2 suffixes to read names (modern convention)
        #-t: Copy RG, BC, and QT tags to FASTQ headers if present
        #--threads 2: Use 2 threads for compression
        #Output File Routing:
        #-1 ~{name}_1.fq.gz: Read 1 of properly paired reads
        #-2 ~{name}_2.fq.gz: Read 2 of properly paired reads
        #-0 >(...): Other reads (neither R1 nor R2, or unpaired)
        #    Reads that don't fit the paired-end model
        #    Could be single-end reads, or reads where pairing information is ambiguous
        #-s >(...): Singleton reads (mate is unmapped/filtered)
        #    Reads whose mate was unmapped or filtered out
        #    One mate passed filters, the other didn't
        #-: Read from stdin (the piped collated BAM)
        
        # Convert BAM to FASTQ
        samtools collate -O -u -n 128 --threads ~{collate_threads} ~{bam} temp \
        | samtools fastq -F 0x900 -n -t --threads 2 \
              -1 ~{name}_1.fq.gz \
              -2 ~{name}_2.fq.gz \
              -0 >(tee >(expr `wc -l` / 4 > count_0.txt) | gzip > ~{name}_0.fq.gz) \
              -s >(tee >(expr `wc -l` / 4 > count_single.txt) | gzip > ~{name}_single.fq.gz) \
              -
    >>>
    
    output {
        File flagstat = '${name}_flagstat.txt'
        File read1 = '${name}_1.fq.gz'
        File read2 = '${name}_2.fq.gz'
        File read0 = '${name}_0.fq.gz'
        File singleton_reads = '${name}_single.fq.gz'
        Int count_0 = read_int('count_0.txt')
        Int count_single = read_int('count_single.txt')
    }
    
    runtime {
        disks: 'local-disk ${runtime_disk_gb} HDD'
        cpu: '${runtime_cpu}'
        memory: '${runtime_memory} GB'
        docker: 'quay.io/biocontainers/samtools:1.17--hd87286a_2'
        preemptible: '${runtime_preemptible}'
    }
}