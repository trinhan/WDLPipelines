## 
## This WDL pipeline runs kallisto pseudo-alignment from single and paired-end fastq files.
## See https://pachterlab.github.io/kallisto/about for more details
##
## 1. Run fastqc (v0.11.8) to check quality of fastqs.
## 2. Run Trimmomatic (v0.39) to remove low-quality basepairs and remove primers/adapter sequence.
## 3. Run kallisto (v0.46.1) to pseudoalign reads to an index/transcriptome. 
##
## Required Inputs:
## - index: Kallisto index for genome of interest. In the form of transcripts.idx
## - fq1: fastq1. Compression state doesn't matter
## - fq2: fastq2. Compression state doesn't matter
## - sample_name
## - adapters for trimmomatic: 
##  String? adapterPath /trimmomatic/adapters/NexteraPE-PE.fa  TruSeq2-PE.fa  TruSeq2-SE.fa  TruSeq3-PE-2.fa  TruSeq3-PE.fa  TruSeq3-SE.fa. Otherwise specify one
## Optional Inputs:
## - preemptiple_attempts: Number of times to allow interruptions on cloud server. Default 3
## - cpu: 1 cpu probably sufficient
## - machine_mem_gb: Memory for each task in GB?  Default: 4 for trimmomatic & fastqc
## - len : estimated fragment length. Required for single-end mode
## - sd : estimated sd for fragment length. Required for single-end mode
##
## Expected Outputs:
## - info: output json file of run commands used
## - abundance: main results of estimated counts and TPMs

version 1.0

workflow KallistoWF {
  input {
    File kallistoindex
    File fq1
    File? fq2
    String sample_name 
    Int cpu = 2
    Int preemptible = 3
    Int mem_gb=8
    ## Inputs for Trimmomatic
    File? adapters
    String? adapterPath
    String? Trimmomaticsettings
    Int iTrimmedReadLen = 36
    ## Inputs for kallisto
    String? len 
    String? sd 
  }

  call fastqc {
    input:
      fastqFile=fq1,
      fastq2=fq2,
      iThreads = cpu
  }


  call RunTrimmomatic { 
    input:
      fq1=fq1,
      fq2=fq2,
      adapters=adapters,
      adapterPath=adapterPath,
      settings=Trimmomaticsettings,
      iTrimmedReadLen=iTrimmedReadLen,
      cpu=cpu
  }

  call quant_kallisto {
    input: 
      index=kallistoindex,
      f1=RunTrimmomatic.fastqOut_R1,
      f2=RunTrimmomatic.fastqOut_R2,
      len=len,
			sd=sd,
      cpu=cpu,
      sample_name=sample_name,
      mem_gb=mem_gb,
      preemptible_attempts=preemptible
  }

  output {
        File infos  = quant_kallisto.info
        File abundance = quant_kallisto.abund
  }
}


task fastqc {
  input {
    File fastqFile
    File? fastq2
    Int iThreads
  }
    
    command {
        ls
        mkdir fastqc_out
        fastqc ~{fastqFile} ~{fastq2} -t ~{iThreads}  -o fastqc_out
    }
    output{
       File zipReport =  "fastqc_out/" + basename(sub(fastqFile,"\\.fastq.gz$","\\_fastqc.zip"))
       File htmlReport = "fastqc_out/" + basename(sub(fastqFile,"\\.fastq.gz$","\\_fastqc.html"))
    }
    runtime{
       docker: "jjkrc/fastqc:0.11.8"
       memory: "4 GB"
       cpu: 2
       disks: "local-disk 200 HDD"
    }
}

task RunTrimmomatic {
  input {
  File fq1
  File? fq2
  Int iTrimmedReadLen
  Int cpu
  File? adapters
  String? adapterPath
  String settings = ":2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15"
  }

  String RunPE = if defined(fq2) then "1" else "0"
  String fq2name = select_first([fq2, ""])
  String adap_def = if defined(adapters) then adapters else select_first([ adapterPath, "/trimmomatic/adapters/TruSeq3-PE.fa"])
  String p_R1 = basename(sub(fq1,".fastq.gz","_trimmed.fastq.gz"))
  String p_R2 = if defined(fq2) then basename(sub(fq2name,".fastq.gz","_trimmed.unpaired.fastq.gz")) else ""
  String u_R1 = basename(sub(fq1,".fastq.gz","_trimmed.unpaired.fastq.gz"))
  String u_R2 = if defined(fq2) then basename(sub(fq2name,".fastq.gz","_trimmed.unpaired.fastq.gz")) else ""
  
  command {

    if [ ~{RunPE} -eq "1" ]; 
    then  
      echo 'run paired send mode'
      java -jar /trimmomatic/trimmomatic.jar PE -threads ~{cpu} -phred33 ~{fq1} ~{fq2} ~{p_R1} ~{u_R1} ~{p_R2} ~{u_R2} ILLUMINACLIP:~{adap_def}~{settings} MINLEN:~{iTrimmedReadLen}
    else
      echo 'run single send mode'
      java -jar /trimmomatic/trimmomatic.jar SE -threads ~{cpu} -phred33 ~{fq1} ~{p_R1} ILLUMINACLIP:~{adap_def}~{settings} MINLEN:~{iTrimmedReadLen}
    fi

  }
  output {
    File fastqOut_R1= "~{p_R1}"
    File? fastqOut_R2= "~{p_R2}"
  }

  runtime {
    docker: "jjkrc/trimmomatic:0.39"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 200 SSD"
  }
}


task RunTrimmomatic_SE {
  input {
  File fq1
  Int iTrimmedReadLen
  Int cpu
  File? adapters
  String? adapterPath
  String settings = ":2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15"
  }

  String adap_def = if defined(adapters) then adapters else select_first([ adapterPath, "/trimmomatic/adapters/TruSeq3-PE.fa"])
  String p_R1 = basename(sub(fq1,".fastq.gz","_trimmed.fastq.gz"))
   
  command {
    java -jar /trimmomatic/trimmomatic.jar SE -phred33 -threads ~{cpu} ~{fq1} ~{p_R1} ILLUMINACLIP:~{adap_def}~{settings} MINLEN:~{iTrimmedReadLen}
  }
  output {
    File fastqOut_R1= "~{p_R1}"
  }

  runtime {
    docker: "jjkrc/trimmomatic:0.39"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 200 SSD"
  }
}

task quant_kallisto {
  input {
    File f1  # File when from input, string when from previous
    File? f2  # File when from input, string when from previous
    File index
    Int cpu
    Int mem_gb
    Int preemptible_attempts
    String sample_name
    String? len = 200
    String? sd = 15
  }

   String SEmode = if defined(f2) then "" else "--single -l ~{len} -s ~{sd} "
   Int disk_size = 3*ceil(size(f1, "GB")+size(f2, "GB")+size(index, "GB"))

  command {

   ## SEcommd=""
   ## if [ ~{SEmode} -eq true ]; then
   ##   SEcommd="--single -l ~{len} -s ~{sd} "
   ## fi
    
    kallisto quant -t ~{cpu} -i ~{index} ~{SEmode} -o out ~{f1} ~{f2}
    mv out/run_info.json ~{sample_name}.json
    sed -re $'s@[|][^\t]*\t@\t@g' out/abundance.tsv > ~{sample_name}.tsv
  
  }

  output {
    File info  = "~{sample_name}.json"
    File abund = "~{sample_name}.tsv"
  }

  runtime {
    docker:"jjkrc/kallisto:0.46.1"
    memory: select_first([mem_gb, 4]) + " GB"
    disk: "local-disk " + disk_size + " HDD"
    cpu: "~{cpu}"
    preemptible: select_first([preemptible_attempts, 3])
  }
}
