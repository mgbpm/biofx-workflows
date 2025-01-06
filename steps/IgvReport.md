# IgvReportFromParsedFASTOutput Task
Takes the parsed FAST output as input and produces a report visualizing the read pile-up for every variant
that is not categorized as "Benign/Likely benign".

## Input Parameters
* File bam_cram - required - the BAM or CRAM file
* File bai_crai - required - the BAM or CRAM index file
* File parsed_fast_output - required - the parsed FAST output, see expected format below 
* String output_basename - optional - the basename that will be used to name the output file, defaults to the parsed FAST output basename with ".igvreport" appended
* File ref_fasta - required - the reference build FASTA file to use for IGV
* File ref_fasta_index - required - the reference build FASTA index file
* Array[File] track_files - optional - a list of files for additional data tracks on the IGV display
* Array[File] track_index_files - optional - corresponding index files for the track files; must have the same basename as
the corresponding track file
* String docker_image - required - the Docker image name and tag
* Int preemptible - optional with a default of 1
* Int colidx_chr - optional with a default of 8. One-based index of the column with Chromosome data.
* Int colidx_start - optional with a default of 9. One-based index of the column with Position Start data.
* Int colidx_end - optional with a default of 10. One-based index of the column with Position End data.

__Parsed FAST Output Format__

The parsed FAST output file is a tab-delimited file with a header row.  The rows following the header can be either:
* Category Header - the category for the subsequent variants (up until the next category header)
* Variant - a variant with values in each of the columns as indicated by the header

There are 3 columns whose position must be provided as an input:
* Chromosome - the chromosome, not including the "chr" prefix
* Start - the variant start position
* End - the variant end position

Additionally, the following columns are expected to be present, but may appear in any position:
* Sym - the gene symbol
* Entrez Gene - the Entrez Gene ID
* Trans - the transcript id
* DNA - the variant DNA change
* AA - the variant amino acid change

## Output Parameters
* File? igv_report - the IGV report HTML file; not present if there are no variants to be included

## Docker Requirements
* Python 3 with modules: igv-reports
