version 1.0

workflow GATKWorflow {

    input {
        File input_cram
        File input_crai
        String sample_id
        String accession_id
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

    String out_path = "outputs/" + sample_id + "/" + test_code + "/supporting/"
    String out_prefix = sample_id + "_" + test_code
    String gvcf = out_path + out_prefix + ".g.vcf"
    String all_calls_vcf = out_path + out_prefix + '.allcalls.vcf'
    String ref_positions_vcf = out_path + out_prefix + '.ref_positions.vcf'
    String all_bases_vcf = out_path + out_prefix + '.allbases.vcf'

    call HaplotypeCallerTask {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            sample_id = sample_id,
            test_code = test_code,
            java_path = java_path,
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

    call GenotypeGVCFsTask {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            sample_id = sample_id,
            test_code = test_code,
            java_path = java_path,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            roi_bed = roi_bed,
            dbsnp = dbsnp,
            dbsnp_vcf_index = dbsnp_vcf_index,
            mgbpmbiofx_docker_image = mgbpmbiofx_docker_image,
            gatk_path = gatk_path 
    }

    output {
        Array[File]+ supporting_files = HaplotypeCallerTask.supporting_files
    }
}


task HaplotypeCallerTask {
    input {
        File input_cram
        File input_crai
        String sample_id
        String test_code
        String java_path
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
        disks: "local-disk 100 SSD"
    }

    output {
        Array[File]+ supporting_files = glob("outputs/" + sample_id + "/" + test_code + "/supporting/*")
    }
}

task GenotypeGVCFsTask {
    input {
        File input_cram
        File input_crai
        String sample_id
        String test_code
        String java_path
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File roi_bed
        File dbsnp
        File dbsnp_vcf_index
        String mgbpmbiofx_docker_image
        String gatk_path
    }


}