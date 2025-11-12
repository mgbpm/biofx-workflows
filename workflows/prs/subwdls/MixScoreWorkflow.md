# PRS Adjustment Workflow

## Input Parameters

| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| Array[File] | input_scores | Yes | PRS scores; each score file's name must contain the PGS ID of the variants weights file used to compute the raw scores (ex: "PGS002236_hmPOS_GRCh38.weights.sscore")  | |
| File | score_weights | Yes | Score weights for each PGS ID | |
| File | output_basename | Yes | Basename for output mix score file | |
| String | ubuntu_docker_image | No | Ubuntu Docker image | "ubuntu:21.10" |

## Output Parameters

| Type | Name | When | Description |
| :--- | :--- | :--- | :--- |
| File | mix_score | Always | Weighted average of input PRS scores |
