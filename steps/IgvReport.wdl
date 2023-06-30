version 1.0

task IgvReportFromParsedFASTOutputTask {
    input {
        File bam_cram
        File bai_crai
        File parsed_fast_output
        String output_basename = sub(basename(parsed_fast_output), "\\.(txt)$", "") + ".igvreport"
        File ref_fasta
        File ref_fasta_index
        String docker_image
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # The bam/cram/ref_fasta index files have to be discoverable by igv-reports
        if [[ "~{bam_cram}" =~ \.bam$ ]]
        then
            [ ! -f "~{bam_cram}.bai" ] && ln -s "~{bai_crai}" "~{bam_cram}.bai"
        elif [[ "~{bam_cram}" =~ \.cram$ ]]
        then
            [ ! -f "~{bam_cram}.crai" ] && ln -s "~{bai_crai}" "~{bam_cram}.crai"
        fi
        [ ! -f "~{ref_fasta}.fai" ] && ln -s "~{ref_fasta_index}" "~{ref_fasta}.fai"

        # Rewrite the parsed fast output into a tab-delimited file 
        #  with the "Benign/Likely benign" rows removed
        awk '
            BEGIN { flag=1; cat=None; FS="\t"; OFS="\t" } 
            {
                if ($1 == "Chr") { print $0, "Category" }
                else if (NF == 1 && $1 == "Benign/Likely benign") { flag=0; cat=$1 } 
                else if (NF == 1 && $1 != "Benign/Likely benign") { flag=1; cat=$1 } 
                else if (NF > 1 && flag == 1) { print $0, cat }
            }
        ' "~{parsed_fast_output}" > igv_sites.txt

        # Check to make certain there are one or more variants in the file
        #   before building the IGV output
        num_sites=$(cat igv_sites.txt | grep -vi ^chr | wc -l)

        if [ $num_sites -gt 0 ]
        then
            create_report igv_sites.txt "~{ref_fasta}" \
                --sequence 1 --begin 2 --end 3 --flanking 50 \
                --info-columns Sym Entrez_Gene Trans DNA AA Category \
                --tracks "~{bam_cram}" \
                --output "~{output_basename}.html"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + ceil(size(bam_cram, "GB") + 15) + " HDD"
        preemptible: preemptible
    }

    output {
        File? igv_report = "~{output_basename}.html"
    }
}