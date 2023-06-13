# QA Testing Workflows
A set of workflows has been created to assist in testing workflow orchestration.
There are five underlying workflow implementations and one or more defined workflows
for each to support testing different configurations of in the same workspace.


## UnitTestWorkflowOne WDL
Accepts 5 string input parameters, which are returned as decorated output parameters and written
to an output file.  Used for the following workflows:
* UnitTestWorkflowSingleParticipant - run workflow with participant as the input entity
* UnitTestWorkflowSingleSample - run workflow with sample as the input entity
* UnitTestWorkflowSingleAlternate - run workflow with a custom entity as the input entity

## UnitTestWorkflowMany WDL
Accepts 5 string array input parameters, which are returned as decorated output parameters and written
to an output file.  Used for the following workflows:
* UnitTestWorkflowMultiParticipant - run workflow with participant set as the input entity
* UnitTestWorkflowMultiSample - run workflow with sample set as the input entity
* UnitTestWorkflowMultiAlternate - run workflow with a custom entity set as the input entity


## E2ETestWorkflowOne
Mocked up single sample CRAM-to-VCF workflow, fetching source files from a remote location.

## E2ETestWorkflowSet
Mocked up sample set CRAM-to-VCF workflow, fetching source files from a remote location.

## TestWorkflowFailure
Workflow that always fails.