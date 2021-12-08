version 1.0

import "https://raw.githubusercontent.com/gatk-workflows/gatk4-germline-snps-indels/master/haplotypecaller-gvcf-gatk4.wdl" as HaplotypeCaller
import "https://raw.githubusercontent.com/trinhan/wgsAlignment/main/pisces_task.wdl" as pisces

workflow runGermlineVariants{
    input {
    File normalBam
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File normalBamIdx
    # a string for the name of the pair under analysis used for naming output files
    String pairName
    # a string for the name of the tumor sample under analysis used for naming output files
    String ctrlName
    # list of read groups to exclude from the analysis in MuTect1 and MuTect_FC tasks
    File refFasta
    # the FASTA file index for the reference genome (see http://www.htslib.org/doc/faidx.html)
    File refFastaIdx
    # the FASTA file dictionary for the reference genome (see https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)
    File refFastaDict
    # an interval list file that contains the locations of the targets
    File targetIntervals
    File? pisces_reference
    Float fracContam 
    String gatk_docker
    }

    String targName=basename(sub(targetIntervals,"\\.interval_list", ""))

    Int normalBam_size  = ceil(size(normalBam,  "G") + size(normalBamIdx,   "G")) 
    Int refFasta_size   = ceil(size(refFasta,   "G") + size(refFastaDict,   "G") + size(refFastaIdx, "G")) 


        # PREPARE FOR SCATTER
    call CallSomaticMutations_Prepare_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            gatk_docker=gatk_docker # takes padded interval file (10bp on each side)
    }

    call CreateFoFN {
        input:
            array_of_files = CallSomaticMutations_Prepare_Task.interval_files,
            fofn_name = targName,
            docker = gatk_docker
    }

    call IntervalToBed {
        input:
            targetIntervals=targetIntervals,
            gatk_docker=gatk_docker,
            output_name = targName
    }

    call Strelka2Germline_Task {
        input: 
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            normalBam_size=normalBam_size,
            refFasta_size=refFasta_size,
            callRegionsBED=IntervalToBed.output_bed,
            callRegionsBEDTBI=IntervalToBed.output_bed_tbi,
            name=pairName
    }

    call pisces.runpisces as runpisces {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaIdx,
            refFastaDict=refFastaDict,
            normalBam=normalBam,
            normalBai=normalBamIdx,
            pairName=pairName,
            pisces_reference=pisces_reference,
            interval=targetIntervals,
            runMode="Germline"   
    }

    call HaplotypeCaller.HaplotypeCallerGvcf_GATK4 as HaplotypeCaller {
        input: 
            input_bam=normalBam,
            input_bam_index=normalBamIdx,
            ref_dict=refFastaDict,
            ref_fasta=refFasta,
            ref_fasta_index=refFastaIdx,
            scattered_calling_intervals_list=CreateFoFN.fofn_list,
            make_gvcf = true,
            make_bamout = false,
            gatk_docker = gatk_docker,
            gatk_path = "/gatk/gatk",
            gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
            samtools_path = "samtools"
    }    

    call Merge_Variants_Germline {
        input:
            ctrlName=ctrlName,
            Haplotype=HaplotypeCaller.output_vcf,
            STRELKA2=Strelka2Germline_Task.strelka2GermlineVCF,
            PISCES_NORMAL=runpisces.normal_variants      
    }

    output {
        # Strelka2Germline
       File strelka2GermlineVCF=Strelka2Germline_Task.strelka2GermlineVCF
       # pisces outputs
       File? pisces_normal_variants=runpisces.normal_variants
       File HaplotypeVcf=HaplotypeCaller.output_vcf
       File HaplotypeVcfTbi=HaplotypeCaller.output_vcf_index
      # merged germline output
       File Merged_germline=Merge_Variants_Germline.MergedGermlineVcf
       File Merged_germlineIdx=Merge_Variants_Germline.MergedGermlineVcfIdx

     }
}

task Strelka2Germline_Task{
    input {
    File normalBam
    File normalBamIdx
    File refFasta
    File refFastaIdx
    File? callRegionsBED
    File? callRegionsBEDTBI
    String name
    String tmpDIR = "strelkaTMP_" + name
    
     # FILE SIZE
    Int normalBam_size
    Int refFasta_size
    String defthreads ="4"

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot ="15"
    String diskGB_buffer ="20"
    String machine_memoryGB ="24"
    String cpu ="1"
}
    # DEFAULT VALUES

    Int command_memoryGB = ceil(machine_memoryGB) - 1
   
    # COMPUTE DISK SIZE
    Int diskGB = ceil(normalBam_size + refFasta_size + diskGB_buffer)


    command {
        mkdir ${tmpDIR} && /usr/local/bin/configureStrelkaGermlineWorkflow.py \
            --bam ${normalBam} \
            --referenceFasta ${refFasta} \
            --exome \
            --runDir ${tmpDIR} --runDir ${tmpDIR} ${"--callRegions " + callRegionsBED} && \
            ${tmpDIR}/runWorkflow.py -m local -j ${defthreads} && \
            mv ${tmpDIR}/results/variants/variants.vcf.gz  ${name}.strelka2.germline.vcf.gz && \
            mv ${tmpDIR}/results/variants/genome.S1.vcf.gz ${name}.strelka2.genome.germline.vcf.gz


    }
    
    runtime{
        docker : "erictdawson/strelka2:2021-Jan-12"
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        cpu            : cpu 
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }
    
    output{
        File strelka2GermlineVCF = "${name}.strelka2.germline.vcf.gz"
    }
}

task CallSomaticMutations_Prepare_Task {
    input {
    # TASK INPUT PARAMS
    File targetIntervals
    File refFasta
    File refFastaIdx
    File refFastaDict

    String nWay = "10"

    # RUNTIME INPUT PARAMS
    String preemptible = "1"
    String diskGB_boot = "15"

    String gatk_docker
}

    parameter_meta {
        nWay : "Number of ways to scatter (MuTect1 and MuTect2)"
        targetIntervals : "a list of genomic intervals over which MuTect1 will operate"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
    }

    command {

        set -euxo pipefail
        mkdir intervalfolder
        gatk SplitIntervals -R ${refFasta} -L ${targetIntervals} --scatter-count ${nWay} -O intervalfolder
        cp intervalfolder/*.interval_list .


    }

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        memory         : "1 GB"
    }

    output {
        Array[File] interval_files=glob("*.interval_list")
    }
}

task CreateFoFN {
  input {
    # Command parameters
    Array[String] array_of_files
    String fofn_name
    # Runtime parameters
    String docker
  }
  command {
    mv ${write_lines(array_of_files)} "${fofn_name}.interval.list"
  }
  output {
    File fofn_list = "${fofn_name}.interval.list"
  }
  runtime {
    docker: docker
    preemptible: 3
  }
}

task IntervalToBed {
    input {
    File targetIntervals
    String? bed_intervallist_extra_args
    String output_name

    # runtime
    String gatk_docker
    String disk_space = "15"
    String preempt = "2"
    String machine_mem = "4"
    }

    Int command_mem = floor(machine_mem) - 1
    String output_name_intervals = output_name + ".bed"

    command {
        set -e
         
        gatk --java-options "-Xmx${command_mem}g" IntervalListToBed \
            -I ${targetIntervals} \
            -O ${output_name_intervals} ${bed_intervallist_extra_args}

        bgzip ${output_name_intervals}

        tabix -s 1 -b 2 -e 3 "${output_name_intervals}.gz"
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preempt
    }

    output {
        File output_bed = "${output_name_intervals}.gz"
        File output_bed_tbi = "${output_name_intervals}.gz.tbi"
    }
}

task Merge_Variants_Germline {

    input {
    # TASK INPUT PARAMS
    File? PISCES_NORMAL
    File? Haplotype
    File? STRELKA2
    String ctrlName
    
    # RUNTIME INPUT PARAMS
    String? preemptible = "2"
    String? diskGB_boot = "10"
    String? diskGB_buffer ="5"
    String? memoryGB ="4"
    String? cpu ="1"
    }

    # DEFAULT VALUES
       
    Int diskGB = ceil(size(Haplotype, "G"))*4  +  diskGB_buffer
    Boolean runS2 =if defined(STRELKA2) then true else false

    command {
        PISCES_pass="~{ctrlName}.Pisces.pass.vcf"
        STRELKA_pass="~{ctrlName}.S2.PASS.vcf.gz"
        HP_unzip="~{ctrlName}.haplo.vcf"
        HP_pass="~{ctrlName}.haplo.pass.vcf"

        MERGED_VCF="~{ctrlName}.P_S2_HP.merged.vcf.gz"
        RENAME_MERGED_VCF="~{ctrlName}.P_S2_HP.mergedGermline.vcf.gz"
 
        ## filter out passed germline variants
        bcftools view -f PASS ${STRELKA2} > $STRELKA_pass
        bcftools view -f PASS ${PISCES_NORMAL} > $PISCES_pass
        ##bgzip $STRELKA_pass
        bgzip $PISCES_pass
        tabix -p vcf $PISCES_pass
        tabix -p vcf $STRELKA_pass.gz

        ## Do the haplotypecaller
        gunzip -c ${Haplotype} > $HP_unzip
        grep -v "0/0"  $HP_unzip > $HP_pass
        bgzip $HP_pass
        tabix -p vcf $HP_pass.gz

        #merge vcfs
        bcftools merge $PISCES_pass.gz $STRELKA_pass.gz $HP_pass.gz -O vcf -o $MERGED_VCF --force-samples
        echo -e "~{ctrlName}.Pisces\n~{ctrlName}.Strelka\n~{ctrlName}.Haplotype\n" > samples.txt
        bcftools reheader -s samples.txt $MERGED_VCF > $RENAME_MERGED_VCF
        tabix -p vcf $RENAME_MERGED_VCF
    }

    runtime {
        docker         : "trinhanne/sambcfhts:v1.13.3"       
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ~{diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
        File MergedGermlineVcf="~{ctrlName}.P_S2_HP.mergedGermline.vcf.gz"
        File MergedGermlineVcfIdx="~{ctrlName}.P_S2_HP.mergedGermline.vcf.gz.tbi"
    }     
}

