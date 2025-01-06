version 1.0

task IgvReportFromParsedFASTOutputTask {
    input {
        File bam_cram
        File bai_crai
        File parsed_fast_output
        String output_basename = sub(basename(parsed_fast_output), "\\.(txt)$", "") + ".igvreport"
        File ref_fasta
        File ref_fasta_index
        Array[File]? track_files
        Array[File]? track_index_files
        String docker_image
        Int preemptible = 1
        Int colidx_chr = 8
        Int colidx_start = 9
        Int colidx_end = 10
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
                if (NR == 1) { print $0, "Parser_Section" } #Output the header row and append the Parser Section column header
                else if (NF == 1 && $1 == "Benign/Likely benign") { flag=0; cat=$1 } #sets section header and set the output flag to off - i.e. ignore Benign/Likely Benign
                else if (NF == 1 && $1 != "Benign/Likely benign") { flag=1; cat=$1 } #sets section header and set the output flag to on - i.e. for all other sections
                else if (NF > 1 && flag == 1) { print $0, cat }
            }
        ' "~{parsed_fast_output}" > igv_sites.txt

        # Check to make certain there are one or more variants in the file
        #   before building the IGV output
        set +e
        num_sites=$(wc -l < igv_sites.txt)
        set -e

        # Prefix chromosome symbols with chr or chrM for MT
        awk -v colidx_chr="~{colidx_chr}" '
            BEGIN {FS=OFS="\t"} 
            NR==1 {print $0; next} 
            {
                if ($colidx_chr == "MT") $colidx_chr = "chrM";
                else $colidx_chr = "chr" $colidx_chr;
            } 
            1' igv_sites.txt > temp && mv temp igv_sites.txt

        # Add spaces after commas so text wrapping in the IGV report table works better
        sed -i -e 's/,\([a-zA-Z0-9]\)/, \1/g' igv_sites.txt

        # Exclude the header from the count
        if [ $num_sites -gt 1 ] 
        then
            track_files=""
            if [ -n "~{sep=' ' track_files}" ]
            then
                mkdir track_files
                mv ~{sep=" " track_files} track_files
                for file in ~{sep=" " track_files}
                do
                    file_name=$(basename $file)
                    track_files="$track_files track_files/$file_name"
                done
            fi
            if [ -n "~{sep=' ' track_index_files}" ]
            then
                mv ~{sep=" " track_index_files} track_files
            fi

            create_report igv_sites.txt "~{ref_fasta}" \
                --sequence ~{colidx_chr} --begin ~{colidx_start} --end ~{colidx_end} --flanking 50 \
                --info-columns Sym "Entrez Gene" Trans DNA AA Parser_Section \
                --tracks "~{bam_cram}" $track_files \
                --output "~{output_basename}.html"
        fi
        exit 0
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

task IgvReportFromGenomePanelsBedTask {
    input {
        File bed_file
        File sample_bam
        File sample_bai
        File ref_fasta
        File ref_fasta_index
        Array[File]? track_files
        Array[File]? track_index_files
        String output_basename 
        String docker_image
        Int preemptible = 1
    }

    command <<<

        set -euxo pipefail

        # The bam/cram/ref_fasta index files have to be discoverable by igv-reports
        if [[ "~{sample_bam}" =~ \.bam$ ]]
        then
            [ ! -f "~{sample_bam}.bai" ] && ln -s "~{sample_bai}" "~{sample_bam}.bai"
        elif [[ "~{sample_bam}" =~ \.cram$ ]]
        then
            [ ! -f "~{sample_bam}.crai" ] && ln -s "~{sample_bai}" "~{sample_bam}.crai"
        fi
        [ ! -f "~{ref_fasta}.fai" ] && ln -s "~{ref_fasta_index}" "~{ref_fasta}.fai"

        # Check to make certain there are one or more variants in the file
        #   before building the IGV output
        num_sites=$(cat "~{bed_file}" | wc -l)

        if [ "$num_sites" -gt 0 ]
        then

            track_files=""
            if [ -n "~{sep=' ' track_files}" ]
            then
                mkdir track_files
                mv ~{sep=" " track_files} track_files
                for file in ~{sep=" " track_files}
                do
                    file_name=$(basename $file)
                    track_files="$track_files track_files/$file_name"
                done
            fi
            if [ -n "~{sep=' ' track_index_files}" ]
            then
                mv ~{sep=" " track_index_files} track_files
            fi

            create_report "~{bed_file}" "~{ref_fasta}" \
                --sequence 1 --begin 2 --end 3 --flanking 50 \
                --info-columns cDNA_Position Predicted_AA \
                --tracks "~{sample_bam}" $track_files \
                --output "~{output_basename}.html"
        fi
    >>>

    runtime {
        docker: "~{docker_image}"
        disks: "local-disk " + ceil(size(sample_bam, "GB") + 15) + " HDD"
        preemptible: preemptible
    }

    output {
        File igv_report = "~{output_basename}.html"
    }

}
