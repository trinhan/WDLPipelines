## STAR fusion pipeline, Jan 2022
## This workflow is an adaptation of the workflow found here: https://github.com/STAR-Fusion/STAR-Fusion
## specifically, using this workflow: https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/Terra-1.10.1/WDL/star_fusion_workflow.wdl
##
## The inputs are:
## - raw fastq files
## 
## The outputs are:
## - fastqc html 
## - Aligned bam file
## - HTSeq counts table
## - STARFusion list of fusions
##
## Key steps:
## 1. Fastqc of files
## 2. Trimmomatic of fastqs
## 3. Run STAR fusion on trimmed fastq
## 4. Index and sort output bam file
## 5. HTSeq on bam file

version 1.0

##import "https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/Terra-1.10.1/WDL/star_fusion_workflow.wdl" as starfusion
import "star_fusion_workflow.wdl" as starfusion
import "FastQC.wdl" as FastQC

workflow STARfusion {

    input {
        String sample_id
        File genome_plug_n_play_tar_gz
        # Inputs required for full pipeline
        File left_fq
        File? right_fq
        String? fusion_inspector
        File gtf #for annotation of genecounts
        Boolean examine_coding_effect
        # Optional inputs for fastqc
        File? adap # input fasta file for adapter sequences, optional
        Boolean runTrimmomatic = true # run trimmomatic on the fastq files
        String? trimmomaticSettings = ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
        # runtime params
        String FusionDocker = "trinityctat/starfusion:1.10.1"
        String samtoolsDocker = "trinhanne/sambcfhts:v1.13.3"
        Int num_cpu = 12
        Float fastq_disk_space_multiplier = 5
        String memory = "50G"
        Float genome_disk_space_multiplier = 5
        Int preemptible = 2
        Float extra_disk_space = 10
        Boolean use_ssd = true
        # small tasks eg. samtools, fastqc etc
        Int smallCPU=1
        Int smallMem=12
    }

    call FastQC.Fastqc as fastqc {
        input:
            f1=left_fq,
            f2=right_fq,
            sampleName=sample_id,
            adap=adap,
            runTrimmomatic=runTrimmomatic,
            settings=trimmomaticSettings,
            cpu=smallCPU
    }

    File f1 = select_first([fastqc.trim_f1, left_fq])
    File f2 = select_first([fastqc.trim_f2, right_fq])

    call starfusion.star_fusion as star_fusion {
        input: 
            left_fq = f1,
            right_fq = f2,
            genome = genome_plug_n_play_tar_gz,
            sample_id = sample_id,
            examine_coding_effect = examine_coding_effect,
            preemptible = preemptible,
            docker = FusionDocker,
            cpu = num_cpu,
            memory = memory,
            extra_disk_space = extra_disk_space,
            fastq_disk_space_multiplier = fastq_disk_space_multiplier,
            genome_disk_space_multiplier = genome_disk_space_multiplier,
            fusion_inspector = fusion_inspector,
            use_ssd = use_ssd
    }

    call samtoolsSortTask {
        input:
            inputBAM = star_fusion.bam,
            sample_id = sample_id,
            dockerImage=samtoolsDocker,
            preemptible=preemptible,
            cpu=smallCPU,
            memoryGB=smallMem
    }

    call HTSeqCount {
        input:
            inputBAM=samtoolsSortTask.outputBAM,
            gtf=gtf,
            preemptible=preemptible,
            sample_id=sample_id,
            cpu=smallCPU,
            memoryGB=smallMem
    }

    output {
        File fusion_predictions = star_fusion.fusion_predictions
        File fusion_predictions_abridged = star_fusion.fusion_predictions_abridged
        File bam = samtoolsSortTask.outputBAM
        File bai = samtoolsSortTask.outputBAI
        File? junction = star_fusion.junction
        File? sj = star_fusion.sj
        File? coding_effect = star_fusion.coding_effect
        Array[File]? extract_fusion_reads = star_fusion.extract_fusion_reads
        File? star_log_final = star_fusion.star_log_final    
        File? fusion_inspector_validate_fusions_abridged = star_fusion.fusion_inspector_validate_fusions_abridged
        File? fusion_inspector_validate_web = star_fusion.fusion_inspector_validate_web
        File? fusion_inspector_inspect_fusions_abridged = star_fusion.fusion_inspector_inspect_fusions_abridged
        File? fusion_inspector_inspect_web = star_fusion.fusion_inspector_inspect_web
        File geneCounts = HTSeqCount.ReadCounts
    }
}

task HTSeqCount {
    input {
        File inputBAM
        File gtf
        String sample_id
        Int preemptible =1
        Int memoryGB = 8
        Int cpu=1
    }
        Int diskSpace=3*ceil(size(inputBAM,"GB")+size(gtf, "GB"))


    command <<<
         htseq-count ~{inputBAM} ~{gtf} -f bam > ~{sample_id}.ReadCounts.txt
    >>>

    runtime {
        docker: "biocontainers/htseq:v0.11.2-1-deb-py3_cv1"
        disks: "local-disk ~{diskSpace} HDD"
        memory: memoryGB + "GB"
        cpu: cpu
        preemptible: preemptible
    }
    output {
        File ReadCounts="~{sample_id}.ReadCounts.txt"
    }
}

task samtoolsSortTask {
    input { 
        File inputBAM
        String sample_id 
        String dockerImage
        Int memoryGB = 16 
        Int cpu = 1
        Int preemptible = 1
    }
        Int diskSpace = 3*ceil(size(inputBAM, "GB"))

    
    command <<<
        samtools sort -o ~{sample_id}.sorted.bam ~{inputBAM} && \
        samtools index ~{sample_id}.sorted.bam
    >>>

    runtime {
        docker: dockerImage
        disks: "local-disk ~{diskSpace} HDD"
        memory: memoryGB + "GB"
        cpu: cpu
        preemptible: preemptible
    }
    output {
        File outputBAM="~{sample_id}.sorted.bam"
        File outputBAI="~{sample_id}.sorted.bam"
    }
}

task ChimFusion {
    input {
      File ChimericJunction
      File genome
      String docker
      String sample_id 
      Int num_cpu
      Float fastq_disk_space_multiplier
      Float genome_disk_space_multiplier
      String? machine_mem_gb
      String? preemptible
      Boolean examine_coding_effect
    }

    Int disk_space_gb = ceil(genome_disk_space_multiplier*size(genome, "GB")+fastq_disk_space_multiplier*size(ChimericJunction, "GB"))

    command <<<

    mkdir -p genome_dir

    tar xf ~{genome} -C genome_dir --strip-components 1

    Input="ChimericInput.Chimeric.out.junction"

    if [[ ~{ChimericJunction} == *.gz ]]; then
        gunzip -c ~{ChimericJunction} > $Input
    else
        cp ~{ChimericJunction} > $Input
    fi

    # Identify the fusions from chimericjunction file

    /usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
             -J $Input \
             --output_dir ~{sample_id} --CPU ~{num_cpu} ~{true='--examine_coding_effect' false='' examine_coding_effect}

    >>>

    output {
      File fusions = "~{sample_id}/star-fusion.fusion_predictions.tsv"
      File fusionsAbridged = "~{sample_id}/star-fusion.fusion_predictions.abridged.tsv"
    }

    runtime {
      docker: docker
      memory: select_first([machine_mem_gb, "50G"])
      cpu: num_cpu
      disks: "local-disk " + disk_space_gb + " HDD"
      preemptible: select_first([preemptible, 3])
    }
}
