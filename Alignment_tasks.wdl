version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.




# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task SamToFastqAndBwaMemAndMba {
  input {
    File input_bam
    String bwa_commandline
    String output_bam_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    Int compression_level
    Int preemptible_tries
    Boolean hard_clip_reads = false
    String gotc_docker  
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(ref_alt, "GiB") + size(ref_amb, "GiB") + size(ref_ann, "GiB") + size(ref_bwt, "GiB") + size(ref_pac, "GiB") + size(ref_sa, "GiB")
  
  # Sometimes the output is larger than the input, or a task can spill to disk.
  # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  
  Float disk_multiplier = 2.5
  
  Int disk_size = ceil(unmapped_bam_size + bwa_ref_size + (disk_multiplier * unmapped_bam_size) + 20)

  command <<<
   BWA_VERSION=$(/usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //')

    set -o pipefail
    set -e

    if [ -z ${BWA_VERSION} ]; then
        exit 1;
    fi
    # set the bash variable needed for the command-line
    bash_ref_fasta=~{ref_fasta}
    # if reference_fasta.ref_alt has data in it,
    if [ -s ~{ref_alt} ]; then
      java -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/~{bwa_commandline} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        REFERENCE_SEQUENCE=~{ref_fasta} \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        ~{true='CLIP_OVERLAPPING_READS=true' false="" hard_clip_reads} \
        ~{true='CLIP_OVERLAPPING_READS_OPERATOR=H' false="" hard_clip_reads} \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="${BWA_VERSION}" \
        PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else reference_fasta.ref_alt is empty or could not be found
   else
      exit 1;
    fi
    
  >>>
  runtime {
    docker: gotc_docker
    preemptible: preemptible_tries
    memory: "14 GiB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

task ConvertToCram {
  input {
    File input_bam
    File ref_fasta
    File ref_fasta_index
    String output_basename
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(2 * size(input_bam, "GiB") + ref_size) + 20

  command <<<
    set -e
    set -o pipefail

    samtools view -C -T ~{ref_fasta} ~{input_bam} | \
    tee ~{output_basename}.cram | \
    md5sum | awk '{print $1}' > ~{output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ~{ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ~{output_basename}.cram
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "~{output_basename}.cram"
    File output_cram_index = "~{output_basename}.cram.crai"
    File output_cram_md5 = "~{output_basename}.cram.md5"
  }
}

# Convert CRAM file to BAM format
task ConvertToBam {
  input {
    File input_cram
    File ref_fasta
    File ref_fasta_index
    String output_basename
  }

  command <<<
    set -e
    set -o pipefail

    samtools view -b -o ~{output_basename}.bam -T ~{ref_fasta} ~{input_cram}

    samtools index ~{output_basename}.bam
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: 3
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk 200 HDD"
  }
  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bam.bai"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  input {
    Array[Float] sizes
    Int preemptible_tries
  }

  command <<<
    python -c "print ~{sep="+" sizes}"
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    preemptible: preemptible_tries
  }
}

task GatherUnsortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/picard/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

workflow SplitLargeReadGroup {

  input {
    File input_bam
    String bwa_commandline
    String output_bam_basename
	
	File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    
    String gotc_docker 
	
    Int compression_level
    Int preemptible_tries
    Int reads_per_file = 50000000
    Boolean hard_clip_reads = false
  }

  call SamSplitter{
    input :
      input_bam = input_bam,
      n_reads = reads_per_file,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    Float current_unmapped_bam_size = size(unmapped_bam, "GiB")
    String current_name = basename(unmapped_bam, ".bam")

    call SamToFastqAndBwaMemAndMba {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = current_name,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries,
        gotc_docker = gotc_docker,
        ref_fasta=ref_fasta,
        ref_fasta_index=ref_fasta_index,
        ref_dict=ref_dict,
        ref_alt=ref_alt,
        ref_amb=ref_amb,
    	ref_ann=ref_ann,
    	ref_bwt=ref_bwt,
    	ref_pac=ref_pac,
    	ref_sa=ref_sa,
        hard_clip_reads = hard_clip_reads
    }

    Float current_mapped_size = size(SamToFastqAndBwaMemAndMba.output_bam, "GiB")
  }

  call SumFloats{
    input:
      sizes = current_mapped_size,
      preemptible_tries = preemptible_tries
  }

  call GatherUnsortedBamFiles {
    input:
      input_bams = SamToFastqAndBwaMemAndMba.output_bam,
      total_input_size = SumFloats.total_size,
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }
  output {
    File aligned_bam = GatherUnsortedBamFiles.output_bam
  }
  meta {
    allowNestedInputs: true
  }
}
