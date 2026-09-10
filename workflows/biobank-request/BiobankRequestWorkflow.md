# Biobank Request Workflow

The Biobank Request Workflow produces a per-request subset of the released MGB
Biobank array, inputed, and/or genomic data.  Given a list of subject IDs and a
list of dataset IDs, it extracts just those subjects from just those datasets
and writes the result to a per-request location in cloud storage, from which it
is delivered to the requesting investigator.

It reads from the *scrubbed* release,
`gs://mgbpm-biobank-data/datasets/current`, so subjects who have withdrawn
consent are already absent from its input.

This workflow is the counterpart of `BiobankScrubWorkflow`, and shares most of
its machinery and its Docker image.  The two differ in the polarity of their
selection: a scrub names the subjects to *exclude* and produces the complement;
a request names the subjects to *include* and produces exactly those.

## Status

**This workflow is a stopgap built on a mistaken premise. It should be
replaced, not completed.**

It was produced under time pressure by repurposing
`BiobankScrubWorkflow`, to meet an IRB requirement that Biobank data
be released as per-request subsets rather than in bulk.  Even though
this approach was able to meet the new requirements with the minimum
of delay, it has proven to be fundamentally wrong.

The resemblance between scrubbing and subsetting is superficial. Both
turn a VCF into a VCF and differ only in which subjects survive, and
that similarity is what made the adaptation quick. But each problem has
features the other lacks.

**The most significant is that a request's subsets may be improper.** A
request can name enough subjects that, for some datasets, every subject
in the dataset has been requested. Those datasets need no subsetting at
all: the deliverable is the released file itself. In the extreme a
request requires no subsetting whatever; more often the delivery is
mixed, part genuine subsets and part verbatim copies of released files.

The scrub never meets this case. There, a file needing no change
requires no action, because it is already in place in `current`. Here, a
dataset needing no subsetting must still be delivered, so "no change"
means "copy", not "do nothing". The two pipelines decide different
questions — *rebuild or leave* against *subset or copy* — and the
deliverable itself has to be conceived differently.

This is not a corner case. Of eleven requests served, three named
64,924, 65,602 and 71,796 subjects, against a cohort of roughly 65,000
unique subjects; for such requests most datasets are wholly contained.

Other differences run the same way. The scrub is triggered by a single
automated notification, runs on a fixed monthly cadence, and treats the
collection uniformly. A request arrives out of a human negotiation, has
no cadence, names an arbitrary set of subjects and datasets, carries an
expiry date, and must reach one named person. None of that is
expressible in machinery designed for the other job.

**A successor should conceptualize this pipeline on its own terms.**
Individual endpoint scripts may well be worth keeping — the subsetting
itself works — but the high-level architecture should not be carried
forward.

This page therefore documents the workflow as an artefact to be
understood, not as a foundation to be built on. It records what the
workflow *does not* do as carefully as what it does. The sections
"Disabled inherited tasks" and "Known gaps" are the substance of it;
read them before "Workflow steps".

It is documented now because its author is leaving, and because until
this commit the WDL existed only as an Agora snapshot and in a single
working copy.

## Input Parameters
| Type | Name | Req'd | Description | Default Value |
| :--- | :--- | :---: | :--- | :--- |
| String | runid | Yes | Identifier for this request; names the output directory under `requestsdir`. By convention the `YYYYMMDD_HHMMSS` timestamp embedded in the Biobank Portal attachment filename | |
| String | operator | Yes | Record only; never read by the workflow. Intended to identify who launched the request. See Known gaps | |
| File | subject_ids | Yes | Text file of Biobank subject IDs, one per line, to be included in the subset | |
| Array[String] | dataset_ids | Yes | The 4-digit dataset IDs to subset, e.g. `["0101", "0102"]` | |
| String | docker_image | Yes | The Docker image supplying the subsetting scripts. Shared with `BiobankScrubWorkflow` | |
| Boolean | keep_shards | No | If true, retain the intermediate per-shard subsets after concatenation instead of purging them | false |
| Int | nbatches | No | Maximum number of batches to partition the subsetting work into | 500 |
| String | requestsdir | No | Cloud storage prefix under which each request's output directory is created | "gs://mgbpm-biobank-data/requests" |
| String | current_datadir | No | Cloud storage prefix of the scrubbed release that supplies the source data | "gs://mgbpm-biobank-data/datasets/current" |
| Int | subset_memory | No | Memory, in GB, for each `SubsetBatch` task | 26 |
| Int | subset_disk_size | No | Local disk, in GB, for each `SubsetBatch` task | 375 |

`docker_image` is declared required, but the WDL carries a commented-out
default naming a `-dev` tag. In the workspace configuration it is supplied as a
literal.

## Output Parameters

**This workflow declares no outputs.** Terra will show an empty OUTPUTS section
for every run.

Results are written directly to cloud storage rather than returned through the
workflow interface. For a given run they appear under:

```
<requestsdir>/<runid>/
```

The individual tasks do declare outputs, which are consumed by later tasks
within the run, but none is promoted to the workflow level.

One consequence deserves emphasis: because there is no `Summarize` task (see
below), **the workflow provides no summary of what it did and no reliable
indication of success**. A green run in Terra is not by itself evidence that
the subset was produced correctly or completely.

## Workflow steps

| # | Task | Scattered | Description |
| :--- | :--- | :---: | :--- |
| 1 | IdentifyPackageCommits | No | Records the git commit of each code package present in the Docker image, for provenance |
| 2 | FindSources | Per `dataset_id` | Determines, for one dataset, which source files under `current_datadir` are needed to satisfy the request |
| 3 | MakeSubsettingBatches | No | Partitions the collected source files into at most `nbatches` batches |
| 4 | SubsetBatch | Per batch | Extracts the requested subjects from each file in the batch, producing subset shards |
| 5 | CollectShards | No | Groups the subset shards by the assembled file each belongs to |
| 6 | ConcatenateShards | Per group | Concatenates each group of shards into a finished VCF; purges the shards unless `keep_shards` is set |

## Disabled inherited tasks

The following tasks are present in the source as commented-out blocks,
inherited from `BiobankScrubWorkflow`. They are listed because their absence is
load-bearing, and because restoring any of them is the most obvious route to
improving this workflow.

| Task | Consequence of its absence |
| :--- | :--- |
| ShowEnvironment | No record of the runtime environment |
| ValidateInputs | Inputs are not checked before work begins |
| MaybeInitializeRundir | No run directory is initialised, so there is no `inputs_history.json` and no record of who launched the run, when, or how often |
| ListDatasetIds | Dataset IDs must be supplied explicitly; they are not discovered |
| CheckBatchedResults | Failures within `SubsetBatch` or `ConcatenateShards` do not gate the steps that follow |
| MakePushReleaseBatches, MakePushShardsBatches, PushScrubbed | No push stage; output remains where `ConcatenateShards` writes it |
| Summarize | No summary document and no aggregate success or failure signal |

## Known gaps

Recorded so that a successor does not have to rediscover them.

- **No success signal.** `Summarize` is disabled, so there is no analogue of
  the scrub's completion check.  How an operator should confirm that a request
  succeeded is undefined.

- **No input validation**, although a purpose-built
  `validate_request_inputs.py` exists in the package. It is called by nothing.

- **No result summarisation**, although `summarize_request_results.py` likewise
  exists and is called by nothing. Both scripts postdate the Docker image in
  use and are absent from it.

- **`operator` is not recorded anywhere.** The workflow never reads it, and
  `MaybeInitializeRundir` — which would have written it to an inputs history —
  is disabled.

- **Chip data is not handled by the deployed image.** `subset_chip.sh` exists
  in the package's working copy but is neither committed nor present in the
  built image. At least one request (dataset `1014`, June 2026) asked for chip
  data only.

- **Version skew.** The Docker image in use predates several commits that were
  written specifically for this workflow. The deployed WDL snapshot is newer
  than the image but does not call the newer scripts.

- **The last steps of the surrounding manual procedure are undocumented** —
  delivery, notification, expiry and deletion. Expiry and deletion matter most,
  since the per-request model exists precisely so that released data can cease
  to be available.

## Related

- `../biobank-scrub/BiobankScrubWorkflow.wdl` — the workflow this one was
  adapted from

- Terra data tables `request` and `message` in the deploying workspace supply
  this workflow's inputs and correlate each request with the email thread that
  originated it
