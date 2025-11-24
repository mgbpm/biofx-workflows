# Run Single PRS Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :--- | :--- | :--- |
| File | query_vcf  | Yes | A gz-compressed VCF file of the samples to be scored | |
| File | adjustment_model_manifest | Yes | JSON file that describes computed model's parameters for either a single weights model | |
| String | condition_code | Yes | Code for condition/disease | |
| Boolean | norename | No | If `true`, do not run `HelperTasks.RenameChromosomes*` tasks | false |
| String | ubuntu_docker_image | No | Ubunutu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | raw_score  | Always | PRS single weights score for the samples in the query VCF |
| File | adjusted_score  | Always | Adjusted PRS single weights score for the samples in the query VCF |
