## Biobank RoR Individual Sample Preparation Workflow
The Biobank Return of Results pipeline is meant to filter Biobank datasets and create VCFs for individual samples within the datasets. This step of the pipeline will prep individual sample VCFs and create a TSV file to upload to the Biobank RoR Terra workspace for the next steps in the pipeline.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[Pair[File,File]] | dataset_files | Yes | An array where each item in the array is a vcf file-index pair; the left being the vcf and the right being the corresponding index file | |
| File | sample_ids_list | Yes | A tsv file containing sample IDs from the dataset | |
| String | dataset | Yes | Name of Biobank dataset | |
| String | dataset_structure | Yes | Whether the dataset is made of joint VCFs or individual sample VCFS; input either "individual" or "joint" | |
| File | target_roi_bed | Yes | BED file containing gene regions; the VCFs in the dataset will be filtered to these regions | |
| Boolean | is_sharded | Yes | Whether or not the input data is sharded | |
| String | bcftools_docker_image | Yes | The name of the bcftools Docker image for VCF annotation | |

## Output Parameters
| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | filtered_vcfs | When filtering input VCFs is applicable | Result of filtering dataset VCFs |
| File | concat_vcf | When concatenating dataset VCFs is applicable; when the dataset is sharded | Result of concatenating dataset VCFs |
| Array[File] | individual_vcfs | When the dataset contains joint VCFs | Result of splitting the dataset VCFs into per sample VCFs |
| File | batch_annotation_input_file | Always | File to use for batch annotation with Alamut and gnomAD; used as input for BiobankRoRAlamutWorkflow |
| File | dataset_sample_table | Always | TSV file containing file paths for individual sample VCFs; upload to Terra to continue the pipeline |