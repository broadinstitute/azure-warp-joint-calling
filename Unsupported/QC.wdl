version 1.0

import "../AzureJointGenotypingTasks.wdl" as Tasks

workflow JointGenotypingQC {
    input {
        File sample_name_map_original_gvcfs
        Array[File] joint_vcf
        Array[File] joint_vcf_index
        Array[File] unpadded_intervals
        String callset_name

        File haplotype_database
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int small_disk
        Int medium_disk
        Int large_disk
        Int huge_disk

        File dbsnp_vcf
        File dbsnp_vcf_index

        # If cross check fingerprints should be scattered, how many gvcfs per shard? Typically set to 1000.
        Int? cross_check_fingerprint_scatter_partition
    }

    Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map_original_gvcfs)
    Int num_gvcfs = length(sample_name_map_lines)
    Array[Array[String]] sample_name_map_lines_t = transpose(sample_name_map_lines)
    Array[String] sample_names_from_map = sample_name_map_lines_t[0]

    call Tasks.CheckSamplesUniqueAndMakeFofn as CheckSamplesUniqueAndMakeFofn {
        input:
            sample_name_map = sample_name_map_original_gvcfs,
            sample_num_threshold = 1
    }

    #TODO: copied from AzureJointGenotypingTasks.wdl, should be moved to a common location

    # CrossCheckFingerprints takes forever on large callsets.
    # We scatter over the input GVCFs to make things faster.
    if (defined(cross_check_fingerprint_scatter_partition)) {
        call Tasks.GetFingerprintingIntervalIndices {
            input:
                unpadded_intervals = unpadded_intervals,
                haplotype_database = haplotype_database
        }

        Array[Int] fingerprinting_indices = GetFingerprintingIntervalIndices.indices_to_fingerprint

        scatter (idx in fingerprinting_indices) {
            File vcfs_to_fingerprint = joint_vcf[idx]
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
            vcf_paths_fofn = write_lines(joint_vcf),
            vcf_index_paths_fofn = write_lines(joint_vcf_index),
            sample_names_from_map_fofn = write_lines(sample_names_from_map),
            gvcf_paths_length = num_gvcfs,
            haplotype_database = haplotype_database,
            output_base_name = callset_name,
            disk = small_disk
        }
    }

    scatter (idx in range(length(joint_vcf))) {
        call Tasks.ValidateVCF {
            input:
                input_vcf = joint_vcf[idx],
                input_vcf_index = joint_vcf_index[idx],
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                calling_interval_list = unpadded_intervals[idx],
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                is_gvcf = false
        }
    }

    output {
        # Output the metrics from crosschecking fingerprints.
        File crosscheck_fingerprint_full = select_first([CrossCheckFingerprintSolo.crosscheck_metrics, GatherFingerprintingMetrics.gathered_metrics])
  }
}