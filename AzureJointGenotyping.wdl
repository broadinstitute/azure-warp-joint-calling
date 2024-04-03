version 1.0

import "AzureJointGenotypingTasks.wdl" as Tasks
import "AzureJointVcfFiltering.wdl" as Filtering


# Joint Genotyping for hg38 Whole Genomes and Exomes (has not been tested on hg19)
workflow JointGenotyping {

  String pipeline_version = "1.6.3"

  input {
    File unpadded_intervals_file
    File targets_interval_list

    String callset_name
    File sample_name_map

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int small_disk
    Int medium_disk
    Int large_disk
    Int huge_disk

    Array[String] snp_recalibration_annotation_values

    File haplotype_database

    File eval_interval_list
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf = dbsnp_vcf
    File dbsnp_resource_vcf_index = dbsnp_vcf_index

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69

    Int? top_level_scatter_count
    Boolean? gather_vcfs
    Int snps_variant_recalibration_threshold = 500000
    Boolean rename_gvcf_samples = true
    Float unbounded_scatter_count_scale_factor = 0.15
    Int gnarly_scatter_count = 10
    Boolean use_gnarly_genotyper = false
    Boolean cross_check_fingerprints = true
    # If cross check fingerprints should be scattered, how many gvcfs per shard? Typically set to 1000.
    Int? cross_check_fingerprint_scatter_partition
  }

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  Array[Array[String]] sample_name_map_lines_t = transpose(sample_name_map_lines)
  Array[String] sample_names_from_map = sample_name_map_lines_t[0]

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  # We allow overriding this default behavior for testing / special requests.
  Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

  Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
  Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?

  call Tasks.CheckSamplesUniqueAndMakeFofn as CheckSamplesUniqueAndMakeFofn {
    input:
      sample_name_map = sample_name_map,
      sample_num_threshold = 5
  }

  call Tasks.SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk,
      sample_names_unique_done = true
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
  	
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call Tasks.ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        # need to provide an example header in order to stream from azure, so use the first gvcf
        header_vcf = CheckSamplesUniqueAndMakeFofn.header_vcf,
        header_vcf_index = CheckSamplesUniqueAndMakeFofn.header_vcf_index,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        disk_size = medium_disk,
        batch_size = 50
    }

    if (use_gnarly_genotyper) {

      call Tasks.SplitIntervalList as GnarlyIntervalScatterDude {
        input:
          interval_list = unpadded_intervals[idx],
          scatter_count = gnarly_scatter_count,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size = small_disk,
          sample_names_unique_done = true
      }

      Array[File] gnarly_intervals = GnarlyIntervalScatterDude.output_intervals

      scatter (gnarly_idx in range(length(gnarly_intervals))) {
        call Tasks.GnarlyGenotyper {
          input:
            workspace_tar = ImportGVCFs.output_genomicsdb,
            interval = gnarly_intervals[gnarly_idx],
            output_vcf_filename = callset_name + "." + idx + "." + gnarly_idx + ".vcf.gz",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp_vcf = dbsnp_vcf
        }
      }

      call Tasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcf_fofn = write_lines(GnarlyGenotyper.output_vcf),
          output_vcf_name = callset_name + "." + idx + ".gnarly.vcf.gz",
          disk_size = large_disk
      }
    }

    if (!use_gnarly_genotyper) {
      call Tasks.GenotypeGVCFs {
        input:
          workspace_tar = ImportGVCFs.output_genomicsdb,
          interval = unpadded_intervals[idx],
          output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          disk_size = medium_disk
      }
    }

    File genotyped_vcf = select_first([TotallyRadicalGatherVcfs.output_vcf, GenotypeGVCFs.output_vcf])
    File genotyped_vcf_index = select_first([TotallyRadicalGatherVcfs.output_vcf_index, GenotypeGVCFs.output_vcf_index])

    call Tasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = genotyped_vcf,
        vcf_index = genotyped_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        targets_interval_list = targets_interval_list,
        disk_size = medium_disk
    }
  }

  #TODO: I suspect having write_lines in the input here is breaking call caching
  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcf_fofn = write_lines(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf),
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size = medium_disk
  }

  call Filtering.JointVcfFiltering as TrainAndApplyVETS {
    input:
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
      input_vcf_idxs = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
      sites_only_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_vcf_idx = SitesOnlyGatherVcf.output_vcf_index,
      targets_interval_list = targets_interval_list,
      annotations = snp_recalibration_annotation_values, #the snp list here is a superset of  the indel list
      resource_vcfs = {"hapmap_resource_vcf": hapmap_resource_vcf,
        "hapmap_resource_vcf_idx": hapmap_resource_vcf_index,
        "omni_resource_vcf": omni_resource_vcf,
        "omni_resource_vcf_idx": omni_resource_vcf_index,
        "one_thousand_genomes_resource_vcf": one_thousand_genomes_resource_vcf,
        "one_thousand_genomes_resource_vcf_idx": one_thousand_genomes_resource_vcf_index,
        "mills_resource_vcf": mills_resource_vcf,
        "mills_resource_vcf_idx": mills_resource_vcf_index},
      output_prefix = callset_name,
      score_runtime_attributes = {"max_retries":1},
      gatk_docker = "mshand/genomesinthecloud:gatk_4_5_0_0"
  }


  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    # For large callsets we need to collect metrics from the shards and gather them later.
    if (!is_small_callset) {
      call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = TrainAndApplyVETS.scored_vcfs[idx],
          input_vcf_index = TrainAndApplyVETS.scored_vcf_idxs[idx],
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size = medium_disk
      }
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
  # HUGE disk was failing in Azure...
  if (is_small_callset) {

    call Tasks.GatherVcfs as FinalGatherVcf {
      input:
        input_vcf_fofn = write_lines(TrainAndApplyVETS.scored_vcfs),
        output_vcf_name = callset_name + ".vcf.gz",
        disk_size = large_disk
    }

    call Tasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size = large_disk
    }
  }

  if (!is_small_callset) {
    # For large callsets we still need to gather the sharded metrics.
    call Tasks.GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = callset_name,
        disk_size = medium_disk
    }
  }

  # CrossCheckFingerprints takes forever on large callsets.
  # We scatter over the input GVCFs to make things faster.
  if (defined(cross_check_fingerprint_scatter_partition)) {
    call Tasks.GetFingerprintingIntervalIndices {
      input:
        unpadded_intervals = unpadded_intervals,
        haplotype_database = haplotype_database,
        disk_size = small_disk
    }

    Array[Int] fingerprinting_indices = GetFingerprintingIntervalIndices.indices_to_fingerprint

    scatter (idx in fingerprinting_indices) {
      File vcfs_to_fingerprint = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx]
    }

    call Tasks.GatherVcfs as GatherFingerprintingVcfs {
      input:
        input_vcf_fofn = write_lines(vcfs_to_fingerprint),
        output_vcf_name = callset_name + ".gathered.fingerprinting.vcf.gz",
        disk_size = medium_disk
    }

    call Tasks.SelectFingerprintSiteVariants {
      input:
        input_vcf = GatherFingerprintingVcfs.output_vcf,
        input_vcf_index = GatherFingerprintingVcfs.output_vcf_index,
        base_output_name = callset_name + ".fingerprinting",
        haplotype_database = haplotype_database,
        disk_size = medium_disk
    }

    # Get partitions by partition number of gvcfs, including any remainder in the last partition
    # Subsetting happens in the CrossCheckFingerprints task
    Array[Int] partitions = range(ceil(num_gvcfs/cross_check_fingerprint_scatter_partition))

    scatter (idx in range(length(partitions))) {
      Int parition_scaled = (partitions[idx] + 1) * cross_check_fingerprint_scatter_partition

      call Tasks.CrossCheckFingerprint as CrossCheckFingerprintsScattered {
        input:
          gvcf_paths_fofn = CheckSamplesUniqueAndMakeFofn.gvcf_paths_fofn,
          gvcf_index_paths_fofn = CheckSamplesUniqueAndMakeFofn.gvcf_index_paths_fofn,
          vcf_paths_fofn = write_lines([SelectFingerprintSiteVariants.output_vcf]),
          vcf_index_paths_fofn = write_lines([SelectFingerprintSiteVariants.output_vcf_index]),
          sample_names_from_map_fofn = write_lines(sample_names_from_map),
          partition_index = parition_scaled,
          partition_ammount = cross_check_fingerprint_scatter_partition,
          gvcf_paths_length = num_gvcfs,
          haplotype_database = haplotype_database,
          output_base_name = callset_name + "." + idx,
          scattered = true,
          disk = small_disk
      }
    }

    call Tasks.GatherPicardMetrics as GatherFingerprintingMetrics {
      input:
        metrics_files = CrossCheckFingerprintsScattered.crosscheck_metrics,
        output_file_name = callset_name + ".fingerprintcheck",
        disk_size = small_disk
    }
  }

  if (!defined(cross_check_fingerprint_scatter_partition)) {

    call Tasks.CrossCheckFingerprint as CrossCheckFingerprintSolo {
      input:
        gvcf_paths_fofn = CheckSamplesUniqueAndMakeFofn.gvcf_paths_fofn,
        gvcf_index_paths_fofn = CheckSamplesUniqueAndMakeFofn.gvcf_index_paths_fofn,
        vcf_paths_fofn = write_lines(TrainAndApplyVETS.scored_vcfs),
        vcf_index_paths_fofn = write_lines(TrainAndApplyVETS.scored_vcf_idxs),
        sample_names_from_map_fofn = write_lines(sample_names_from_map),
        gvcf_paths_length = num_gvcfs,
        haplotype_database = haplotype_database,
        output_base_name = callset_name,
        disk = small_disk
    }
  }

  # Get the metrics from either code path
  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  # Get the VCFs from either code path
  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else TrainAndApplyVETS.scored_vcfs
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else TrainAndApplyVETS.scored_vcf_idxs

  output {
    # Metrics from either the small or large callset
    File detail_metrics_file = output_detail_metrics_file
    File summary_metrics_file = output_summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)

    # Output the interval list generated/used by this run workflow.
    Array[File] output_intervals = SplitIntervalList.output_intervals

    # Output the metrics from crosschecking fingerprints.
    File crosscheck_fingerprint_check = select_first([CrossCheckFingerprintSolo.crosscheck_metrics, GatherFingerprintingMetrics.gathered_metrics])
  }
  meta {
    allowNestedInputs: true
  }
}
