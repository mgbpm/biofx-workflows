# Run PRSmix Workflow

This workflow will calculate a PRSmix score for samples within a query VCF. There are two options for methods to calculate the PRSmix score:

1) Score the query VCF with each variant weights file, mix the raw scores, and then adjusting the mix score with a single PRSmix model
2) Score the query VCF with each variant weights file, adjust each raw score with its own model, and then mix the adjusted scores

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | query_vcf  | Yes | A gz-compressed VCF file of the samples to be scored | |
| Array[File] | adjustment_model_manifests | Yes | Array of JSON files that describe computed model's parameters for either single weights models (to mix PRS scores after adjustment) or one PRSmix model (to mix PRS scores before adjustment) | |
| String | condition_code | Yes | Code for condition/disease | |
| Boolean | mix_before_adjustment | Yes | If `true`, mix raw scores before adjusting them |  |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | final_score  | Always | Adjusted PRS single weights score or mix score for the samples in the query VCF |
