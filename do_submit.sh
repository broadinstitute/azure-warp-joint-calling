zip /tmp/tasks.zip AzureJointGenotypingTasks.wdl
echo "{}" > /tmp/empty.json
cromshell submit AzureJointGenotyping.wdl AzureJointGenotyping.22samples.inputs.json /tmp/empty.json /tmp/tasks.zip
