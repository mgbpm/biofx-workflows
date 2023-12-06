# LMMGAP Workflow
The LMMGAP Workflow processes GTC files from an Illumina genotype array batch.  It converts
the GTC files to VCF files and performs a series of QC checks on the resulting data.

## Input Parameters
| Type          | Name                          | Req'd | Description                                                                   | Default Value                                                                   |
|:--------------|:------------------------------|:-----:|:------------------------------------------------------------------------------|:--------------------------------------------------------------------------------|
| String        | gcp_project_id                | No    | The GCP project id to fetch connection secrets from                           | "mgb-lmm-gcp-infrast-1651079146"                                                |
| String        | workspace_name                | Yes   | The name of the workspace the workflow is running in                          |                                                                                 |
| String        | orchutils_docker_image        | No    | The orchutils Docker image to use                                             | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/orchutils:20230828"           |
| String        | lmmgap_docker_image           | No    | The lmmgap Docker image to use                                                | "gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/lmmgap:20231206"              |
| String        | batch_id                      | No    | A prefixed value for naming ouptput and intermediate files                    | Random GUID                                                                     |
| String        | project_file_location         | Yes   | A path to a file that supplies details pertaining to the samples for analysis |                                                                                 |
| Array[String] | subject_ids                   | Yes   | A list of subject ids of the samples to be processed, indexed to sample_ids   |                                                                                 |
| Array[String] | sample_ids                    | Yes   | A list of identifiers of the samples to be processed                          |                                                                                 |
| Array[String] | sample_genders                | Yes   | A list of sample genders, indexed to sample_ids                               |                                                                                 |
| Array[String] | sample_data_locations         | Yes   | A list of sample data location paths, indexed to sample_ids                   |                                                                                 |
| File          | genome_reference_fa           | Yes   | The reference assembly FASTA file                                             |                                                                                 |
| File          | genome_reference_fai          | Yes   | The reference assembly FASTA index file                                       |                                                                                 |
| File          | genome_reference_dict         | Yes   | The reference assembly DICT file                                              |                                                                                 |
| String        | genome_reference_name         | No    | The reference assembly name                                                   | "hg19"                                                                          |
| File          | reference_sample_project_file | No    | The reference sample project detail CSV file                                  | "gs://lmm-reference-data/lmmgap/27202778065-8_ReferenceProjectDetailReport.csv" |
| String        | reference_sample_gender       | No    | The gender of the reference sample                                            | "2"                                                                             |
| String        | reference_subject_id          | No    | The subject id of the sample used as a reference                              | "NA12878-A10"                                                                   |
| String        | reference_sample_id           | No    | The identifier of the sample used as a reference                              | "204339030100-R01C01"                                                           |  
| File          | reference_sample_gtc_file     | No    | Reference sample GTC file                                                     | "gs://lmm-reference-data/lmmgap/NA12878-A10_204339030100-R01C01.gtc"            |
| Float         | call_rate_cut_off             | No    | Call rate for QC determination                                                | 0.99                                                                            |
| File          | bpm_manifest_file             | No    | Path to bpm manifest file                                                     | "gs://lmm-reference-data/lmmgap/GDA-8v1-0_A5.bpm"                               |
| File          | csv_manifest_file             | No    | Path to csv manifest file                                                     | "gs://lmm-reference-data/lmmgap/GDA-8v1-0_A5.csv"                               |
| File          | excluded_site_file            | No    | Path to excluded probes file                                                  | "gs://lmm-reference-data/lmmgap/GDA-8v1-0_A5.probes-to-exclude.tsv"             |
| Boolean       | save_working_files            | No    | Whether to include working files as workflow outputs                          | true                                                                            |

## Output Parameters
| Type        | Name                     | Description                                                |
| :-----------|:-------------------------|:-----------------------------------------------------------|
| File        | final_report             | The final report |
| File        | master_file              |  |
| Array[File] | sample_vcfs              | One VCF file for each sample |
| Array[File] | sample_vcf_idxs          | An index file for each VCF file |
| File        | merged_vcf               | The merged VCF file containing all samples |
| File        | merged_bed               | The merged BED file |
| File        | merged_bed_bim           |  |
| File        | merged_bed_map           |  |
| File        | merged_bed_hh            |  |
| File        | merged_bed_ped           |  |
| File        | merged_bed_fam           |  |
| File        | split_bed                |  |
| File        | split_bed_hh             |  |
| File        | split_bed_bim            |  |
| File        | split_bed_fam            |  |
| File        | genome_qc                | PLINK genome QC output |
| File        | genome_qc_hh             |  |
| File        | sexcheck_qc              | PLINK sex-check QC output |
| File        | sexcheck_qc_hh           |  |
| File        | missing_qc_imiss         | PLINK missing QC output |
| File        | missing_qc_lmiss         |  |
| File        | missing_qc_hh            |  |
| Array[File] | working_gtc_to_vcf_files |  |
| Array[File] | working_fixhdr_vcf_files |  |
| Array[File] | working_sorted_vcf_files |  |
| Array[File] | working_fixhet_vcf_files |  |
| Array[File] | working_log_files        |  |