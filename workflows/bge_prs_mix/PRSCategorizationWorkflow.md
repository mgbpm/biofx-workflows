# PRS Categorization Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | prs_scores | Yes | PRS scores for each condition that is desired in the final summary | |
| File | sample_ids | Yes | All IDs for samples replated to the PRS scores | |
| String | ubuntu_docker_image | No | Ubuntu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| Array[File] | individuals_risk_summaries | Always | Risk summary report with risk percentile and bins for each sample |