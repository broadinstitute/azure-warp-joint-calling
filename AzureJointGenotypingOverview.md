# JointGenotyping on Azure

| Pipeline | Date Updated | Documentation Authors | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [AzureJointGenotyping](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotyping.wdl) | February, 2024 | Kaylee Mathews & Megan Shand | Please file issues in GitHub |

## Introduction to the AzureJointGenotyping workflow

The [AzureJointGenotyping workflow](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotyping.wdl) is an open-source, cloud-optimized pipeline that implements joint variant calling and filtering using using GATK and Microsoft Azure. The pipeline calls the [Variant Extract-Train-Score (VETS) subworkflow](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointVcfFiltering.wdl) to score variant annotations.

The AzureJointGenotyping pipeline can be configured to run using one of the following GATK joint genotyping methods:

* **[GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/21905118377755)** (default method) performs joint genotyping on GVCF files stored in GenomicsDB and pre-called with HaplotypeCaller.
* **[GnarlyGenotyper](https://gatk.broadinstitute.org/hc/en-us/articles/21904951112091)** performs scalable, “quick and dirty” joint genotyping on a set of GVCF files stored in GenomicsDB and pre-called with HaplotypeCaller.

The pipeline takes in a sample map file listing GVCF files produced by HaplotypeCaller or DRAGEN version 3.7.8 in GVCF mode and creates a filtered VCF file (with index) containing genotypes for all samples present in the input VCF files. All sites that are present in the input VCF file are retained. Filtered sites are annotated as such in the FILTER field. If you are new to VCF files, see the [file type specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

> **Note**
> The pipeline is adapted from the [WARP JointGenotyping workflow](https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping/README), but is not subject to the same testing requirements as WARP pipelines. 


## Set-up

The AzureJointGenotyping pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform.

### Inputs

The [AzureJointGenotyping workflow](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotyping.wdl) is adapted from the [WARP JointGenotyping workflow](https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping/README) and requires all the same inputs with the following exceptions:

| Parameter name | Description | Type |
| --- | --- | --- |
| targets_interval_list | Required input. | File |
| snp_recalibration_tranche_values | Removed from pipeline. | Array[String] |
| indel_recalibration_tranche_values | Removed from pipeline. | Array[String] |
| indel_recalibration_annotation_values | Removed from pipeline. | Array[String] |
| vqsr_snp_filter_level | Removed from pipeline. | Float |
| vqsr_indel_filter_level | Removed from pipeline. | Float |
| snp_vqsr_downsampleFactor | Removed from pipeline. | Int |
| use_allele_specific_annotations | Removed from pipeline. | Boolean |
| run_vets | Removed from pipeline. | Boolean |
| scatter_cross_check_fingerprints | Removed from pipeline; scattering during fingerprinting is determined by `cross_check_fingerprint_scatter_partition`. | Boolean |
| cross_check_fingerprint_scatter_partition | Optional integer specifying the number of samples to include in each partition for scattering during fingerprinting; recommended value is “1000”; fingerprinting will be performed without scattering if no value is passed to the pipeline. | Int |
| sample_name_map | Path to file containing the sample names (first column; example: “NA12878”) and the Azure cloud path of the individual GVCF files (second column; example: “az://sc-74cc28aa-fa7c-4712-8b3e-7eb784790bec@lzb25a77f5eadb0fa72a2ae7.blob.core.windows.net/path_to_file.NA12878.vcf.gz”). | String | 


## Azure JointGenotyping tasks and tools

Overall, the AzureJointGenotyping workflow:

1. Splits the input interval list and imports GVCF files.
2. Performs joint genotyping using GATK GenotypeGVCFs (default) or GnarlyGenotyper.
3. Creates single site-specific VCF and index files.
4. Creates and applies a variant filtering model using VETS.
5. Collects variant calling metrics.
6. Checks fingerprints.

The [AzureJointGenotyping workflow](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotyping.wdl) imports individual “tasks,” also written in WDL script. 

The [AzureJointGenotyping workflow](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotyping.wdl) is adapted from the [WARP JointGenotyping workflow](https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping/README) and calls all the same tasks with the following exceptions:


| Task | Tool | Software | Description | 
| --- | --- | --- | --- | 
| [CheckSamplesUniqueAndMakeFofn as CheckSamplesUniqueAndMakeFofn](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | bash | bash | Renamed from `CheckSamplesUnique`; checks that there are more than 50 unique samples in `sample_name_map` and generates necessary sample map files for Azure. |
| [JointVcfFiltering as TrainAndApplyVETS](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointVcfFiltering.wdl) | ExtractVariantAnnotations, TrainVariantAnnotationsModel, ScoreVariantAnnotations | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Default method for variant filtering; calls the `JointVcfFiltering.wdl` subworkflow to extract variant-level annotations, trains a model for variant scoring, and scores variants. | 
| [IndelsVariantRecalibrator](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Removed from the pipeline. | 
| [SNPsVariantRecalibratorCreateModel](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Removed from the pipeline. | 
| [SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Removed from the pipeline. | 
| [GatherTranches as SNPGatherTranches](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | GatherTranches | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Removed from the pipeline. | 
| [SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Removed from the pipeline. | 
| [ApplyRecalibration](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | ApplyVQSR | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Removed from the pipeline. |
| [GetFingerprintingIntervalIndices](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | IntervalListTools | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprint_scatter_partition` is defined, gets and sorts indices for fingerprint intervals; otherwise the task is skipped. |
| [GatherVcfs as GatherFingerprintingVcfs](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | GatherVcfsCloud | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprint_scatter_partition` is defined, compiles the fingerprint VCF files; otherwise the task is skipped. |
| [SelectFingerprintSiteVariants](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | SelectVariants | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprint_scatter_partition` is defined, selects variants from the fingerprint VCF file; otherwise the task is skipped. |
| [PartitionSampleNameMap](https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl) | bash | bash | Removed from the pipeline. |
| [CrossCheckFingerprint as CrossCheckFingerprintsScattered](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | CrosscheckFingerprints | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprint_scatter_partition` is defined, checks fingerprints for the VCFs in the scattered partitions and produces a metrics file; otherwise the task is skipped. |
| [GatherPicardMetrics as GatherFingerprintingMetrics](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | bash | bash | If `cross_check_fingerprint_scatter_partition` is defined, combines the fingerprint metrics files into a single metrics file; otherwise the task is skipped. |
| [CrossCheckFingerprint as CrossCheckFingerprintSolo](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotypingTasks.wdl) | CrosscheckFingerprints | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprint_scatter_partition` is not defined, checks fingerprints for the single VCF file and produces a metrics file; otherwise the task is skipped. |


## Outputs

The [AzureJointGenotyping workflow](https://github.com/broadinstitute/azure-warp-joint-calling/blob/main/AzureJointGenotyping.wdl) is adapted from the [WARP JointGenotyping workflow](https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping/README) and outputs all the same files.