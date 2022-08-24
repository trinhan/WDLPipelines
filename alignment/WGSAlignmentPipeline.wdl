version 1.0

## This is a copy of the Broad Institute WGS Alignment pipeline, modified 12/11/2021
##
## Changes compared to the original:
## - Split uBAM if too large, and parallelise mapping
## - Change preemptible options if the output file is too big 
## 
## Copyright Broad Institute, 2019
## 
## This WDL pipeline implements data pre-processing according to the GATK Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Default Options:
## - docker: broadinstitute/gatk: 4.2.6.1
## - gotc: broadinstitute/genomes-in-the-cloud:2.3.1-1512499786
## - cutoff_size: if raw fastq is greater than [10GB], it will be split
## - reads_per_file: divide the fastq into files with [50 000 000] reads each. For 30x, this is ~20 shards
## - preempt_cutoff: if the file is greater than [110 GB] preemptibles will not be used for sorting and marking dups
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.
##
## TO INCLUDE:
## 1. Fastq to uBAM? (Pipeline ready)
## 2. QC METRICS?
## 3. BAM TO CRAM?

import "Alignment_tasks.wdl" as Alignment
import "AggregatedBamQCmod.wdl" as QC
import "BamCramConvert.forGitHub.wdl" as Bam2Cram

# WORKFLOW DEFINITION 
workflow PreProcessingForVariantDiscovery_GATK4 {
  input {
    String sample_name
    String ref_name

    File flowcell_unmapped_bams_list
    String unmapped_bam_suffix
  
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    File wgs_coverage_interval_list

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
    Int compression_level = 5
  
    String gatk_docker = "broadinstitute/gatk:4.2.6.1"
    String gatk_path = "/gatk/gatk"
    String gotc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    String gotc_path = "/usr/gitc/"
    String python_docker = "python:2.7"  
    Boolean hard_clip_reads = false
    
    Int preemptible_tries = 3
    Float cutoff_size = 10
    Int read_length = 150
    Int reads_per_file = 50000000
    Int preempt_cutoff = 110

    Boolean CheckPreVal = false
    Boolean compressCram
  }
    String base_file_name = sample_name + "." + ref_name
    String RGchecksum_md5 = base_file_name + ".RGchecksum.md5"

    Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)
    

 # Get the version of BWA to include in the PG record in the header of the BAM produced 
 # by MergeBamAlignment. 
  call GetBwaVersion {
    input: 
      docker_image = gotc_docker,
      bwa_path = gotc_path,
      preemptible_tries = preemptible_tries
  }

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in flowcell_unmapped_bams) {

   Float unmapped_bam_size=size(unmapped_bam, "GiB")
    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(unmapped_bam, unmapped_bam_suffix)
  
  if (unmapped_bam_size > cutoff_size) {
      # Split bam into multiple smaller bams,
      # map reads to reference and recombine into one bam
      call Alignment.SplitLargeReadGroup as SplitLargeReadGroup {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          bwa_version=GetBwaVersion.version,
          output_bam_basename = bam_basename + ".aligned.unsorted",
          compression_level = compression_level,
          ref_fasta=ref_fasta,
          ref_fasta_index=ref_fasta_index,
          ref_dict=ref_dict,
          gotc_docker=gotc_docker,
          gotc_path=gotc_path,
          bwa_path=gotc_path,
          reads_per_file=reads_per_file,
          gatk_docker=gatk_docker,
          gatk_path=gatk_path,
          preemptible_tries=preemptible_tries,
          ref_alt=ref_alt,
          ref_amb=ref_amb,
          ref_ann=ref_ann,
          ref_bwt=ref_bwt,
          ref_pac=ref_pac,
          ref_sa=ref_sa
      }
    }

if (unmapped_bam_size <= cutoff_size) {
    # Map reads to reference
    call Alignment.SamToFastqAndBwaMem as SamToFastqSingle {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = bam_basename,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries,
        gotc_docker = gotc_docker,
        gotc_path=gotc_path,
        bwa_path=gotc_path,
        ref_fasta=ref_fasta,
        ref_fasta_index=ref_fasta_index,
        ref_dict=ref_dict,
        ref_alt=ref_alt,
        ref_amb=ref_amb,
        ref_ann=ref_ann,
        ref_bwt=ref_bwt,
        ref_pac=ref_pac,
        ref_sa=ref_sa
     }
     # Mege the bam files togeher
     call Alignment.MergeBamAlignment as MergeBam {
      input:
        unmapped_bam = unmapped_bam,
        bwa_commandline = bwa_commandline, # put in a call to figure out bwa_version
        bwa_version=GetBwaVersion.version,
        aligned_bam=SamToFastqSingle.output_bam,
        output_bam_basename=bam_basename,
        ref_fasta=ref_fasta,
        ref_fasta_index=ref_fasta_index,
        ref_dict=ref_dict,
        compression_level=4,
        preemptible_tries=preemptible_tries,
        gatk_docker=gatk_docker,
        gatk_path=gatk_path
     }
     
     }
   File output_aligned_bam = select_first([MergeBam.output_bam, SplitLargeReadGroup.aligned_bam])
   Float mapped_bam_size = size(output_aligned_bam, "GiB")  

    }
    
    call Alignment.SumFloats as SumFloats {
    input:
      sizes = mapped_bam_size,
      preemptible_tries = 3
  }
    
    Boolean data_too_bigs = SumFloats.total_size > preempt_cutoff

    call MarkDuplicates {
    input:
      input_bams = output_aligned_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",   
      docker_image = gatk_docker,
      total_input_size=SumFloats.total_size,
      gatk_path = gatk_path,
      compression_level = compression_level,
      preemptible_tries = if data_too_bigs then 0 else preemptible_tries
  }

  # Sort aggregated+deduped BAM file and fix tags
  ## agg_large_disk
  call SortAndFixTags {
    input:
      input_bam = MarkDuplicates.markdup_output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      preemptible_tries = if data_too_bigs then 0 else preemptible_tries,
      compression_level = compression_level
  }
  
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      docker_image = python_docker,
      preemptible_tries = preemptible_tries
  }
  
   Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
   Int potential_bqsr_divisor = num_of_bqsr_scatters - 25
   Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1
   
   scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortAndFixTags.sort_output_bam,
        input_bam_index = SortAndFixTags.sort_output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        bqsr_divisor=bqsr_divisor,
        preemptible_tries = preemptible_tries
    }  
  }
  
  
   call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      preemptible_tries = preemptible_tries
  }   
  
  
   Int num_of_bqsr_scatters2 = length(CreateSequenceGroupingTSV.sequence_grouping_with_unmapped)
   Int potential_bqsr_divisor2 = num_of_bqsr_scatters2 - 25
   Int bqsr_divisor2 = if potential_bqsr_divisor2 > 1 then potential_bqsr_divisor2 else 1
 
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    
    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.sort_output_bam,
        input_bam_index = SortAndFixTags.sort_output_bam_index,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        bqsr_divisor=bqsr_divisor2,
        preemptible_tries = preemptible_tries
    }
  }
  
  
  Float agg_bam_size = size(ApplyBQSR.recalibrated_bam, "GiB")
  
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      total_input_size=agg_bam_size,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }


  ## Add a QC step here:

  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = base_file_name + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list, ## need this
      read_length = read_length,
      preemptible_tries = preemptible_tries
  }

  call QC.AggregatedBamQC as AggregatedBamQC { 
    input:
      base_recalibrated_bam = GatherBamFiles.output_bam,
      base_recalibrated_bam_index = GatherBamFiles.output_bam_index,
      base_name = base_file_name,
      BSQR_md5 = RGchecksum_md5,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
  }

if (compressCram){
  call Bam2Cram.BamToCram as BamToCram {
    input:
      input_file = GatherBamFiles.output_bam,
      base_file_name = base_file_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      duplication_metrics = MarkDuplicates.duplicate_metrics,
      chimerism_metrics = AggregatedBamQC.agg_alignment_summary_metrics,
      ToCram = true,
      agg_preemptible_tries=0,
      CheckPreVal=CheckPreVal
  }
}

  # Outputs that will be retained when execution is complete  
  output {
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File analysis_ready_bam = GatherBamFiles.output_bam
    File analysis_ready_bam_index = GatherBamFiles.output_bam_index
    File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5

    File wgs_metrics = CollectWgsMetrics.metrics

    File read_group_alignment_summary_metrics = AggregatedBamQC.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = AggregatedBamQC.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = AggregatedBamQC.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = AggregatedBamQC.read_group_gc_bias_summary_metrics
    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5
    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = AggregatedBamQC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = AggregatedBamQC.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File? cram_out=BamToCram.output_file
    File? crai_out=BamToCram.output_index
    File? cram_md5=BamToCram.output_md5_file
    File? cram_validate_report=BamToCram.validate_report
  } 
   
   
}


##TASKS
# Get version of BWA
task GetBwaVersion {
  input {
    Float mem_size_gb = 1
    Int preemptible_tries
    String docker_image
    String bwa_path
  }  

  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    BWA_VERSION=~{bwa_path}bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
  }
  output {
    String version = read_string(stdout())
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  input {
    File input_bam
    String output_bam_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
  
    Int compression_level
    Int preemptible_tries
    Float mem_size_gb = 32

    String docker_image
    String gatk_path
  }
    
    Int command_mem_gb_sort = ceil(mem_size_gb) - 1
    Int command_mem_gb_fix = ceil((mem_size_gb - 1)/10)
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
 
  command {
    set -o pipefail

    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_sort}G" \
      SortSam \
      --INPUT ~{input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_fix}G" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ~{ref_fasta}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File sort_output_bam = "~{output_bam_basename}.bam"
    File sort_output_bam_index = "~{output_bam_basename}.bai"
    File sort_output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
  
    Int compression_level
    Int preemptible_tries
    Float mem_size_gb = 7.5
    Int additional_disk = 20
    Float total_input_size
    String docker_image
    String gatk_path
  }
    Int command_mem_gb = ceil(mem_size_gb) - 2
    Float md_disk_multiplier = 3
    Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk

    
    
 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
      MarkDuplicates \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{output_bam_basename}.bam \
      --METRICS_FILE ~{metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb}  GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File markdup_output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task CreateSequenceGroupingTSV {
 input {
    File ref_dict  
  
    Int preemptible_tries
    Float mem_size_gb = 2

    String docker_image
  }
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
  
    python <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
  
    Int preemptible_tries
    Float mem_size_gb = 6
    Int bqsr_divisor
    
    String docker_image
    String gatk_path
  }
  
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Float dbsnp_size = size(dbSNP_vcf, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB")/bqsr_divisor)*1.2 + ref_size + dbsnp_size) + 20
  
  Int command_mem_gb = ceil(mem_size_gb) - 2

  command { 
    ~{gatk_path} --java-options "-Xms~{command_mem_gb}G" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbSNP_vcf} \
      --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename

    Int compression_level
    Int preemptible_tries
    
    Float mem_size_gb = 3
    Float total_input_size

    String docker_image
    String gatk_path
  }
  
  Int command_mem_gb = ceil(mem_size_gb) - 1
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
      GatherBamFiles \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int preemptible_tries
    Float mem_size_gb = 4
    Int additional_disk = 20
    String docker_image
    String gatk_path
    
    Int bqsr_divisor
  }
  Int command_mem_gb = ceil(mem_size_gb) - 1
  
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") * 3/bqsr_divisor ) + ref_size) + additional_disk

  Int memory_size = "4"


  command {  
    ~{gatk_path} --java-options "-Xms~{command_mem_gb}G" \
      ApplyBQSR \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{output_bam_basename}.bam \
      -L ~{sep=" -L " sequence_group_interval} \
      -bqsr ~{recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
  }
}


task GatherBqsrReports {
 input { 
   Array[File] input_bqsr_reports
   String output_report_filename

   Int preemptible_tries
   Int disk_size = 15
   Float mem_size_gb = 4

   String docker_image
   String gatk_path
  }
  Int command_mem_gb = ceil(mem_size_gb) - 1

  command {
    ~{gatk_path} --java-options "-Xms~{command_mem_gb}G" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bqsr_report = "~{output_report_filename}"
  }
}


