## Strelka2 workflow
## Requirements
## Bam: tumor if run in tumor mode, normal Required
## runMode: "Tumour" or "Germline"
## runGerm: Do you want to run Germline (if in tumor mode?)
## strelka_config file: See
## bed file (.gz & .gz.tbi) of regions to profile
## test removal of config file


version 1.0

workflow Strelka2FullWorkflow {
    input {
    File refFasta
    File refFastaIdx
    File? tumorBam
    File? tumorBamIdx
    File normalBam
    File normalBamIdx
    File? callRegionsBED
    File? callRegionsBEDTBI
    Int? threads
    String pairName
    File? strelka_config
    String runMode ="Tumour" ## or "Germline"
    Boolean runGerm = true
    }

    Int tumorBam_size=ceil(size(tumorBamIdx, "GB")+ size(tumorBam, "GB"))
    Int normalBam_size=ceil(size(normalBamIdx, "GB")+ size(normalBam, "GB"))
    Int refFasta_size=ceil(size(refFastaIdx, "GB")+ size(refFasta, "GB"))

    if (runMode=="Tumour"){
        call Strelka2SomaticTask {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            callRegionsBED=callRegionsBED,
            callRegionsBEDTBI=callRegionsBEDTBI,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            refFasta_size=refFasta_size,
            name=pairName,
            config=strelka_config,
            defthreads=threads
        }
    }

    if (runMode=="Germline" || runGerm){
        call Strelka2GermlineTask {
        input: 
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            normalBam_size=normalBam_size,
            refFasta_size=refFasta_size,
            callRegionsBED=callRegionsBED,
            callRegionsBEDTBI=callRegionsBEDTBI,
            normalBam_size=normalBam_size,
            refFasta_size=refFasta_size,
            name=pairName,
            threads=threads
        }
    }

    output {
        File? strelkaGermlineVCF=Strelka2GermlineTask.strelkaGermlineVCF
        File? strelka2SomaticSNVs = Strelka2SomaticTask.strelka2SomaticSNVs
        File? strelka2SomaticIndels = Strelka2SomaticTask.strelka2SomaticIndels
    }

}



task Strelka2SomaticTask {
    input {
    # TASK INPUT PARAMS
    File? tumorBam
    File? tumorBamIdx
    File normalBam
    File normalBamIdx
    File refFasta
    File refFastaIdx
    File? config
    File? callRegionsBED
    File? callRegionsBEDTBI
    String name
   
    String tmpDIR = "strelkaTMP_" + name

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size
    Int? defthreads ="4"

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot = "15"
    String diskGB_buffer ="20"
    String machine_memoryGB ="25"
    String cpu ="1"
    }
    # DEFAULT VALUES
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + diskGB_buffer)

    command {

        set -euxo pipefail
        
        mkdir ${tmpDIR} && /usr/local/bin/configureStrelkaSomaticWorkflow.py \
            --normalBam ${normalBam} \
            --tumorBam ${tumorBam} \
            --referenceFasta ${refFasta} \
            --exome \
            --runDir ${tmpDIR} ${"--callRegions " + callRegionsBED} && \
            ${tmpDIR}/runWorkflow.py -m local -j ${defthreads} && \
            mv ${tmpDIR}/results/variants/somatic.snvs.vcf.gz ${name}.strelka.somatic.snvs.vcf.gz && \
            mv ${tmpDIR}/results/variants/somatic.indels.vcf.gz ${name}.strelka.somatic.indels.vcf.gz
    }

    runtime {
        docker         : "erictdawson/strelka2:2021-Jan-12"
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File strelka2SomaticSNVs = "${name}.strelka2.somatic.snvs.vcf.gz"
        File strelka2SomaticIndels = "${name}.strelka2.somatic.indels.vcf.gz"
    }
}

task Strelka2GermlineTask{
    input {
    File refFasta
    File refFastaIdx
    File normalBam
    File normalBamIdx
    File? callRegionsBED
    File? callRegionsBEDTBI
    Int? threads ="4"
    String name
    String tmpDIR = "strelkaTMP_" + name
    Int normalBam_size
    Int refFasta_size
    
    String preemptible ="1"
    String diskGB_boot = "15"
    String diskGB_buffer ="20"
    String machine_memoryGB ="25"
    String cpu ="1"
    }
    
   
    # COMPUTE DISK SIZE

    Int diskGB = ceil(normalBam_size + refFasta_size + diskGB_buffer)

 


    command {
        mkdir ${tmpDIR} && /usr/local/bin/configureStrelkaGermlineWorkflow.py \
            --bam ${normalBam} \
            --referenceFasta ${refFasta} \
            --exome ${"--callRegions " + callRegionsBED} && \
            --runDir ${tmpDIR} && \
            ${tmpDIR}/runWorkflow.py -m local -j ${threads} && \
            mv ${tmpDIR}/results/variants/variants.vcf.gz ${name}.strelka.germline.vcf.gz 
    }
    
    runtime {
        docker : "erictdawson/strelka2:2021-Jan-12"
         bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }
    
    output { 
        File strelkaGermlineVCF = "${name}.strelka.germline.vcf.gz"
    }
}
