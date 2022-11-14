# azure-warp-joint-calling
Exploration of Joint Calling Exomes with Cromwell on Azure using WARP Joint Calling Pipeline

The Joint Calling pipeline and tasks WDLs are originally from the WARP repository, and modified here for running on Azure at small scale (~10k exomes).

The example input files are for joint calling of 22 exomes that are part of the AnVIL 1000 Genomes public data set.

## Description of Files

Pipeline:
 - `AzureJointGenotyping.wdl` - Main WDL
 - `AzureJointGenotypingTasks.wdl` - Subtasks included from main WDL
 
Test Data / Scipts:
 - `AzureJointGenotyping.22samples.inputs.json` - Pipeline inputs for 22-sample exome data set
 - `azure_sample_map.22samples.txt` - sample map for 22-sample exome data set
 - `do_submit.sh` - example script to zip up task and use Cromshell to submit workflow
 
 ## How To Run
 
The prototype was run against a CromwellOnAzure set up under subscription `8201542-gvs-dev` in the resource group `test-cromwell-on-azure-rg`.  Storage is all under the storage account `coa6fb978543c0ccf` where all the associated Azure storage containers are located.

In order to enable the use of Cromshell, one can ssh into the Cromwell server with ssh port forwarding back to the local machine on port 8000 with

```
ssh -L 8000:localhost:8000 vmadmin@coa-130787859d14ef.eastus.cloudapp.azure.com
```

And then configure cromshell to talk to localhost:8000 for Cromwell.

All the pipeline inputs (gVCFs, reference files, and the sample map input) must be in the `inputs` storage container and referenced.

Finally the `do_submit.sh` shell script zips up the tasks and uses cromshell to submit the pipeline.
 
 
