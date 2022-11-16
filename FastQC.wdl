## script to run fastqc and 

version 1.0

workflow Fastqc {
  input {
    File f1
    File? f2
    String sampleName
    Int? cpu = 1
    # Adaptor (check documentation)
    File? adap
    Boolean runTrimmomatic
    String settings = ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
  }
    String mode = if defined(f2) then "PE" else "SE"

      call fastqc as fastq1 {
        input: f1=f1,
               f2=f2,
               cpu=cpu,
               mode=mode
      }

  if (runTrimmomatic){
        call trim_PE {
            input: f1=f1,
                   f2=f2,
                   cpu=cpu,
                   sampleName=sampleName,
                   adap=adap,
                   settings=settings,
                   mode=mode
        }

      call fastqc as fastq2 {
        input: f1=trim_PE.out_p1,
               f2=trim_PE.out_p2,
               cpu=cpu,
               mode=mode
      }
    }

  

  output {
    File pre_html = fastq1.html 
    File pre_zip = fastq1.zip
    File? post_html = fastq2.html 
    File? post_zip = fastq2.zip
    File? trim_f1 = trim_PE.out_p1
    File? trim_f2 = trim_PE.out_p2
  }
}

task trim_PE {
  input {
    File f1
    File? f2
    Int? cpu=1
    Int? machine_mem_gb
    Int? requested_disk
    Int? preemptible_attempts
    File? adap
    String sampleName
    String mode
    String settings = ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
  }
  Int disk_space_gb=ceil(3*(size(f1, "GB")+size(f2, "GB")))
  File adap_def = select_first([adap,"/opt/trimmomatic/adapters/TruSeq3-PE.fa"])
  String p1 = sampleName + "Trim_R1.fq.gz"
  String p2 = sampleName + "Trim_R2.fq.gz"
  String u1 = sampleName + "Trim_unpaired_R1.fq.gz"
  String u2 = sampleName + "Trim_unpaired_R2.fq.gz"
  command <<<
      if [ ~{mode} == "PE" ]; then
      trimmomatic PE -threads ~{cpu} -phred33 ~{f1} ~{f2} \
            ~{p1} ~{u1} ~{p2} ~{u2} ILLUMINACLIP:~{adap_def}~{settings}
      else
      trimmomatic SE -threads ~{cpu} -phred33 ~{f1} ~{p1} ILLUMINACLIP:~{adap_def}~{settings}
      fi
  >>>

  output {
    File out_p1 = "~{p1}"
    File? out_p2 = "~{p2}"
    File out_u1 = "~{u1}" 
    File? out_u2 = "~{u2}" 
  }

  runtime {
    docker: "bmwlee/usyd_trimmomatic:0.36"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    cpu: cpu
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
}


task fastqc {
  input {
    File f1
    File? f2
    Int? cpu
    Int? machine_mem_gb
    Int? preemptible_attempts
    String mode
    Int? requested_disk
  }
  # Disk space calculation
  Float f2size = if (mode == "PE") then size(f2, "GB") else 0.0
  Int predicted_disk = ceil(3*(size(f1, "GB")+f2size))
  Int disk_space_gb = select_first([requested_disk,predicted_disk])

  command <<<

    f2comm=""
    if [ mode == "PE" ]; then
    f2comm="~{f2}"
    fi

    fastqc -t ~{cpu} --outdir $PWD ~{f1} $f2comm
  >>>

  output {
    File html = glob("*html")[0]
    File zip  = glob("*zip")[0]
  }

  runtime {
    docker: "biocontainers/fastqc:v0.11.5_cv2"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    cpu: cpu
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
}
