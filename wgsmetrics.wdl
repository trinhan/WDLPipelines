## CollectWgsmetrics By interval or WGS
## 
## This workflow runs wgs metrics based on list of interval files submitted.
## This text file must contain the paths for all intervals that are run in parallel
##
## Required:
## - bam and bam index
## - sampleName
## - reference genome and index
## - Intervals: txt file containing list of regions to assess in scatter mode. 
##   Default: gs://kccg-cb-annotation-files/BedIntervals/hg38/wgs_contigs/wgs_1_22_chrom_list
## 
## Optional:
## - ByChromosome: whether to run on each chr individually (true) or whole genome (false)
## - FastMode: turn on fast mode for wgsmetrics? can be done for wgs analysi
## - GATK docker image

version 1.0
workflow WgsMetrics {
  input {
    File bam
    File bam_index
    String sampleName
    File refFasta
    File refFastaIdx
    File Intervals
    Boolean ByChromosome = false
    Boolean FastMode = false
    String docker
    String opts = ""
  }

if (ByChromosome == true){

  Array[File] intervals_list = read_lines(Intervals)
  scatter (chromList in intervals_list){
    String chr = basename(sub(chromList, "\\.interval_list", ""))
    call runWgsMetrics as runWgsMetricsChrom {
       input:
           bam=bam,
           bam_index=bam_index,
           sampleName=sampleName,
           refFasta=refFasta,
           refFastaIdx=refFastaIdx,
           Intervals=chromList,
           FastMode=FastMode,
           docker=docker,
           region=chr,
           opts=opts
    }
  }

    call compressFiles {
      input:
        outputsChrom = runWgsMetricsChrom.stats,
        sampleName = sampleName,
        docker=docker
    }
}

if (ByChromosome == false ){
 call runWgsMetrics as runWgsMetricsAll {
       input:
           bam=bam,
           bam_index=bam_index,
           sampleName=sampleName,
           refFasta=refFasta,
           refFastaIdx=refFastaIdx,
           FastMode=FastMode,
           docker=docker,
           region="wgs",
           opts=opts
    }
}

  File Output = select_first([compressFiles.outputsTar, runWgsMetricsAll.stats])

output {
  File outputs = Output
}

}


task compressFiles {
  input{
    Array[File] outputsChrom
    String sampleName
    String docker
  }
  command <<<

  tar czf ~{sampleName}.wgs.by.chrom.tar.gz ~{outputsChrom}
  >>>

output {
  File outputsTar = "~{sampleName}.wgs.by.chrom.tar.gz"
}
  runtime {
    docker: docker
    preemptible: "3"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }

}

task runWgsMetrics {
  input{
    File bam
    File bam_index
    String sampleName
    File refFasta
    File refFastaIdx
    File? Intervals
    Boolean FastMode = false
    String docker
    String region
    String opts
  }
    Int diskspace = 2*ceil(size(bam, "GB")+size(bam_index, "GB")+size(refFasta, "GB")+size(refFastaIdx, "GB"))

  command <<<

  gatk CollectWgsMetrics -I ~{bam} -O ~{sampleName}.~{region}.txt -R ~{refFasta} \
   --USE_FAST_ALGORITHM ~{FastMode} ~{"--INTERVALS " + Intervals} ~{opts}

  >>>
  output {
   File stats="~{sampleName}.~{region}.txt"
  }
  runtime {
    docker: docker
    preemptible: "3"
    memory: "5 GB"
    disks: "local-disk ~{diskspace} HDD"
  }
}