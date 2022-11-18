# Mass General Brigham Personalized Medicine Bioinformatic Workflows
This repository contains all of the WDL-based bioinformatic workflows and tasks created by 
Mass General Brigham Personalized Medicine for processing genomic data.  It also contains 
utility and orchestration tasks.  The reusable steps and workflows are organized into a 
directory structure, outlined below by way of example:

```
steps/                                - Contains all re-usable workflow steps, which may be WDL tasks and workflows
  ...wdls                             - WDL files, with one workflow and/or many tasks per file
  ...mds                              - Workflow and task documentation files, one for each WDL file
workflows/                            - Contains one sub-directory for each MGBPM workflow
  wgs_b38/                            - Contains all artifacts for the workflow
    WGSSingleSampleBuild38.wdl        - The workflow definition
    WGSSingleSampleBuild38.md         - Workflow documentation, especially inputs and outputs
  prs_b38/
    ...similar to wgs...
.dockstore.yml                        - File that defines all the workflows to be imported by Dockstore
```

# WDL Authoring Guidelines
_ALL CONTENTS OF THIS REPOSITORY WILL BE PUBLISHED ON GITHUB_

_CREDENTIALS OR OTHER CONFIDENTIAL INFORMATION ARE STRICTLY FORBIDDEN_

The following guidelines should be followed for WDL development:
* All tasks that execute MGBPM authored code must accept container image name as an input parameter that defaults to `gcr.io/mgb-lmm-gcp-infrast-1651079146/mgbpmbiofx/base:latest`
* The WDLs should be designed to be consumed as black box functions; i.e. all input parameters and data are defined as inputs and all generated data of interest to consumers is defined as an output
* Workflows should be reference build specific; reusable steps should be reference build agnostic
* Workflows should always require reference data (other than the reference build), sample data and external system connection information to be passed as parameters

When workflow WDLs are ready for production use, the relevant commit will be tagged with workflow specific tags to avoid collisions.  For example, `WGSSingleSampleBuild38_2.5.0`.  Versioning for workflow WDLs should adhere to the principles of [semantic versioning](https://semver.org/).

The master copy of this repository is hosted in Azure DevOps.  All new content is mirrored to GitHub and workflows subsequently published in Dockstore.
