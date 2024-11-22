# PRS Adjustment Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| File | condition_name | Yes | Name of condition/disease | |
| Array[File] | raw_scores | Yes | Raw PRS scores | |
| File | score_weights | Yes | Score weights for each PGS ID | |
| String | ubuntu_docker_image | No | Ubuntu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | prs_mix_raw_score | Always | Per sample weighted average of raw PRS scores |
| File | sample_ids | Always | Sample IDs corresponding to the samples in the PRS score file |
