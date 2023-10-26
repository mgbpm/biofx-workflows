version 1.0

workflow GATKWorkflow {

    input {
        File input_cram
        File input_crai
        String sample_id
        String test_code
        String java_path = "/usr/lib/jvm/java-8-openjdk-amd64/bin/java"
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        String gatk_path = "/gatk/gatk"
        String mgbpmbiofx_docker_image
    }

    String out_path = "outputs/"
    String out_prefix = sample_id + "_" + test_code
    String gvcf = out_path + out_prefix + ".g.vcf"
    String all_calls_vcf = out_path + out_prefix + '.allcalls.vcf'
    String ref_positions_vcf = out_path + out_prefix + '.ref_positions.vcf'
    String all_bases_vcf = out_path + out_prefix + '.allbases.vcf'

    call HaplotypeCallerTask {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gatk_path = gatk_path,
            out_path = out_path,
            gvcf = gvcf,
            all_calls_vcf = all_calls_vcf
    }

    call CreateRefSitesVCFTask {
        input:
            gvcf_file = HaplotypeCallerTask.gvcf_file,
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            ref_positions_vcf = ref_positions_vcf,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            out_path = out_path
    }

    call SortVCFTask {
        input:
            gatk_path = gatk_path,
            java_path = java_path,
            out_path = out_path,
            ref_positions_vcf_file = CreateRefSitesVCFTask.ref_positions_vcf_file,
            all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file,
            all_bases_vcf = all_bases_vcf,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
    }

    output {
        File all_calls_vcf_file = HaplotypeCallerTask.all_calls_vcf_file
        File all_calls_vcf_idx_file = HaplotypeCallerTask.all_calls_vcf_idx_file
        File gvcf_file = HaplotypeCallerTask.gvcf_file
        File gvcf_idx_file = HaplotypeCallerTask.gvcf_idx_file
        File ref_positions_vcf_file = CreateRefSitesVCFTask.ref_positions_vcf_file
        File all_bases_vcf_file = SortVCFTask.all_bases_vcf_file
        File all_bases_vcf_idx_file = SortVCFTask.all_bases_vcf_idx_file
    }
}


task HaplotypeCallerTask {
    input {
        File input_cram
        File input_crai
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        String mgbpmbiofx_docker_image
        String gatk_path
        String out_path
        String gvcf
        String all_calls_vcf
        Int disk_size = ceil(size(input_cram, "GB") + size(dbsnp, "GB")) + 10
    }
    
    command <<<
        set -euxo pipefail

        mkdir -p ~{out_path}

        ~{gatk_path} --java-options "-Xmx20G" \
        HaplotypeCaller \
        --input ~{input_cram} \
        --output ~{gvcf} \
        --reference ~{reference_fasta} \
        --intervals ~{roi_bed} \
        --dbsnp ~{dbsnp} \
        -mbq 10 \
        -stand-call-conf 30 \
        -A QualByDepth \
        -A FisherStrand \
        -A RMSMappingQuality \
        -A MappingQualityZero \
        -A StrandOddsRatio \
        -A DepthPerAlleleBySample \
        --read-filter MappingQualityReadFilter \
        --minimum-mapping-quality '17' \
        --read-filter MappingQualityNotZeroReadFilter \
        -ERC BP_RESOLUTION
        
        ~{gatk_path} --java-options "-Xmx20G" \
        GenotypeGVCFs \
        --variant ~{gvcf} \
        --output ~{all_calls_vcf} \
        --reference ~{reference_fasta} \
        --intervals ~{roi_bed} \
        --dbsnp ~{dbsnp} \
        -A QualByDepth \
        -A FisherStrand \
        -A RMSMappingQuality \
        -A MappingQualityZero \
        -A StrandOddsRatio \
        -A DepthPerAlleleBySample \
        --read-filter MappingQualityReadFilter \
        --minimum-mapping-quality '17' \
        --read-filter MappingQualityNotZeroReadFilter
        
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "24 GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        bootDiskSizeGb: 15
    }

    output {
        File all_calls_vcf_file = "~{all_calls_vcf}"
        File all_calls_vcf_idx_file = "~{all_calls_vcf}.idx"
        File gvcf_file = "~{gvcf}"
        File gvcf_idx_file = "~{gvcf}.idx"
    }
}

task CreateRefSitesVCFTask {
    input{
        File gvcf_file
        File all_calls_vcf_file
        String ref_positions_vcf
        String mgbpmbiofx_docker_image
        String out_path
        Int disk_size = ceil(size(all_calls_vcf_file, "GB") + size(gvcf_file, "GB")) + 10
    }

    command <<<
        set -euxo pipefail

        mkdir -p ~{out_path}

        $MGBPMBIOFXPATH/biofx-*/bin/create_ref_sites_vcf.py \
        -g "~{gvcf_file}" \
        -c "~{all_calls_vcf_file}" \
        -o "~{ref_positions_vcf}"
        
    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        bootDiskSizeGb: 15
    }

    output {
        File ref_positions_vcf_file = "~{ref_positions_vcf}"
    } 
}

task SortVCFTask {
    input {
        String gatk_path
        String java_path
        File ref_positions_vcf_file
        File all_calls_vcf_file 
        String all_bases_vcf
        String mgbpmbiofx_docker_image
        String out_path
        Int disk_size = ceil(size(all_calls_vcf_file, "GB") + size(ref_positions_vcf_file, "GB")) + 10
    }

    command <<<
        set -euxo pipefail

        mkdir -p ~{out_path} 
        
        ~{java_path} -Xms12g -Xmx40g \
        -jar ~{gatk_path}.jar SortVcf \
        -I ~{ref_positions_vcf_file} \
        -I ~{all_calls_vcf_file} \
        -O ~{all_bases_vcf}

    >>>

    runtime {
        docker: "~{mgbpmbiofx_docker_image}"
        memory: "44 GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        bootDiskSizeGb: 15
    }   

    output {
        File all_bases_vcf_file = "~{all_bases_vcf}"
        File all_bases_vcf_idx_file = "~{all_bases_vcf}.idx"
    }
}