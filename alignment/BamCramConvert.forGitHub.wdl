version 1.0

# Local Import
#import "Utilities.wdl" as Utils
#import "Qc.wdl" as QC

# Git URL Import

import "https://raw.githubusercontent.com/trinhan/wgsAlignment/main/Alignment_tasks.wdl?token=ABVSYKGLUML35DJTE76BS73BAYBLW" as Alignment

workflow BamToCram {

  input {
    File input_file
    #File input_cram
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? duplication_metrics
    File? chimerism_metrics
    String base_file_name
    Int agg_preemptible_tries
    Boolean ToCram = true
    Boolean CheckPreVal = false
  }


  # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
  # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
  Float max_duplication_in_reasonable_sample = 0.30
  Float max_chimerism_in_reasonable_sample = 0.15

  # Convert the final merged recalibrated BAM file to CRAM format
  if (ToCram) {
    call Alignment.ConvertToCram as ConvertToCram {
      input:
        input_bam = input_file,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        output_basename = base_file_name,
        preemptible_tries = agg_preemptible_tries
    }
  }

  if (ToCram==false) {
    call Alignment.ConvertToBam as ConvertToBam {
      input:
        input_cram = input_file,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        output_basename = base_file_name,
        preemptible_tries = agg_preemptible_tries      
    }
  }

  File outputfile=select_first([ConvertToCram.output_cram, ConvertToBam.output_bam])
  File outputindex=select_first([ConvertToCram.output_cram_index, ConvertToBam.output_bam_index])
  String validatesuffix=if ToCram then "cram" else "bam"
  File output_md5=select_first([ConvertToCram.output_cram_md5, ConvertToBam.md5_bam])

  if (CheckPreVal && defined(duplication_metrics) && defined(chimerism_metrics)){
    call CheckPreValidation {
      input:
        duplication_metrics = duplication_metrics,
        chimerism_metrics = chimerism_metrics,
        max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
        max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample,
        preemptible_tries = agg_preemptible_tries
    }
  }

  Boolean? is_outlier_data=select_first([CheckPreValidation.is_outlier_data, false])

  Array[String] ignorethis = if ToCram then ["MISSING_TAG_NM"] else ["null"]
      

  call ValidateSam {
    input:
      input_bam = outputfile,
      input_bam_index = outputindex,
      report_filename = base_file_name + validatesuffix+ ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ignore =ignorethis,
      max_output = 1000000000,
      is_outlier_data=is_outlier_data,
      preemptible_tries = agg_preemptible_tries
  }


  ## check the QC metrics:
  output {
     File output_file= outputfile
     File output_index = outputindex
     File output_md5_file = output_md5
     File validate_report = ValidateSam.report
  }
}

#####################
## PUT TASKS HERE!!##
#####################


task ValidateSam {
  input {
    File input_bam
    File? input_bam_index
    String report_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int? max_output
    Array[String]? ignore
    Boolean? is_outlier_data
    Int preemptible_tries
    Int memory_multiplier = 1
    Int additional_disk = 20
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + additional_disk

  Int memory_size = ceil(7 * memory_multiplier)
  Int java_memory_size = (memory_size - 1) * 1000

  command {
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File report = "~{report_filename}"
  }
}

task CheckPreValidation {
  input {
    File? duplication_metrics
    File? chimerism_metrics
    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    Int preemptible_tries
  }

  command <<<
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python <<CODE

    import csv
    with open('duplication.csv') as dupfile:
      reader = csv.DictReader(dupfile, delimiter='\t')
      for row in reader:
        with open("duplication_value.txt","w") as file:
          file.write(row['PERCENT_DUPLICATION'])
          file.close()

    with open('chimerism.csv') as chimfile:
      reader = csv.DictReader(chimfile, delimiter='\t')
      for row in reader:
        with open("chimerism_value.txt","w") as file:
          file.write(row['PCT_CHIMERAS'])
          file.close()

    CODE

  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    preemptible: preemptible_tries
    memory: "2 GiB"
  }
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  }
}
