### This pipeline takes a list of variants and performs:
## 1. blat realignment and removes failed variants (hg19 only, requires genome_bit file)
## 2. perform indel realignment using abra2
## 3. Re-calls variants using strelka2, pisces, mutect1, mutect2 after indel realignment
## 4. Outputs merged variants file which has been filtered by Panel of Normals and 
## Required inputs:
## bam files, reference files, 


version 1.0

import "abra2.wdl" as abra2
import "SNV-MultiCaller.wdl" as SomaticVC
import "blat_runner.wdl" as blat

workflow variantJudgement {
    input {
    # Required variant list
    File merged_maf
    File vcf
    File vcfidx
    File indeltargets
    # Bam files
    File tumorBam
    File tumorBamIdx
    File? normalBam
    File? normalBamIdx
    String pairName
    String? ctrlName
    String caseName
    # Ref Genome
    File refFasta
    File refFastaIdx
    File refFastaDict
    File targetIntervals
    File indeltargets
    # To run indelrealignment: true or false
    Boolean runAbra2Realign
    Boolean runBlatRealign
    # Annotation Files
    File? genome_bit
    File gnomad
    File gnomadidx
    File PoN
    File PoNidx    
    File? variants_for_contamination
    File? variants_for_contamination_idx
    File? DB_SNP_VCF
    File? DB_SNP_VCF_IDX
    File? cosmicVCF
    File knownSVs
    File knownSVsidx

    Float fracContam
    String runMode
    String gatk_docker
    Int tumorBam_size
    File strelka_config
    String refGenome
    Int minCallerSupport

    Boolean callM2
    Boolean callStrelka
    Boolean callVardict
    Boolean callM1
    Boolean callPisces = false

    }

    Boolean runBlat = if ( defined(genome_bit) && refGenome =="hg19" && runBlatRealign == true) then true else false
    File gene_bit = select_first([genome_bit, "NULL"])

    if (runBlat){
    call blat.blat as blat {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            MAF=merged_maf, 
            pairName=pairName,
            genome_bit=gene_bit,
            tumorBam_size=tumorBam_size,
            refGenome=refGenome
    }
    }

    if (runAbra2Realign) {
    call abra2.runabra2 as runabra2 {
        input:
            normalBam=normalBam,
            normalIdx=normalBamIdx,
            tumorBam=tumorBam,
            tumorIdx=tumorBamIdx,
            refFasta=refFasta,
            vcf=knownSVs,
            vcfidx=knownSVsidx,
            sample_name=pairName,
            refFastaIdx=refFastaIdx,
            targets=indeltargets
    }

    call SomaticVC.runVariantCallers as somaticVC {
        input:
            normalBam=runabra2.Abra2Normal,
            normalBamIdx=runabra2.Abra2NormalIdx,
            tumorBam=runabra2.Abra2Tumour,
            tumorBamIdx=runabra2.Abra2TumourIdx,
            refFasta=refFasta,
            pairName=pairName,
            caseName=caseName,
            ctrlName=ctrlName,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            gnomad=gnomad,
            gnomad_idx=gnomadidx,
            m2_pon=PoN,
            m2_pon_idx=PoNidx,
            variants_for_contamination=variants_for_contamination,
            variants_for_contamination_idx=variants_for_contamination_idx,
            DB_SNP_VCF=DB_SNP_VCF,
            DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
            cosmicVCF=cosmicVCF, # can this be updated?
            fracContam=fracContam,
            runMode=runMode,
            gatk_docker=gatk_docker,
            minCallerSupport=minCallerSupport,
            strelka_config=strelka_config,
            callM2=callM2,
            callStrelka=callStrelka,
            callPisces=callPisces,
            callVardict=callVardict,
            callM1=callM1
    }
    }

    call PassJudgement {
        input:
            blat_reject=blat.rejectMaf,
            abra2_calls=somaticVC.Combined_raw_variants_gz,
            vcf=vcf,
            m2_pon=PoN,
            gnomad=gnomad,
            pairName=pairName,
            runMode=runMode,
            refGenome=refGenome
    }

     output {
        File? blat_maf=blat.passMaf
        File? blat_reject=blat.rejectMaf
        File? abra2_tum_bam=runabra2.Abra2Tumour
        File? abra2_tum_bam_idx=runabra2.Abra2TumourIdx
        File? abra2_norm_bam=runabra2.Abra2Normal
        File? abra2_norm_bam_idx=runabra2.Abra2NormalIdx
        File? merged_abra2_calls=somaticVC.Combined_raw_variants_gz
        File judgement_vcf=PassJudgement.pass_judgement
        File judgement_maf=PassJudgement.pass_maf
    }
}


task PassJudgement {
    input {
        File? blat_reject
        File? abra2_calls
        File vcf
        File? gnomad
        File? m2_pon
        String pairName
        String refGenome
        String runMode

        # RUNTIME INPUT PARAMS
        String? preemptible = "2"
        String? diskGB_boot = "10"
        String? diskGB_buffer ="5"
        String? memoryGB ="10"
        String? cpu ="1"
    }

    # DEFAULT VALUES
       
    Int diskGB = ceil(size(m2_pon, "G")+size(gnomad, "G")+size(vcf, "G"))*4 +  diskGB_buffer

    command {
        # Use an R script for this?
        # change blat to vcf file
        # bcftools isec for PoN
        
        Rscript /app/MergeVCFRecall.R --vcffile ${vcf} --outputstr ${pairName} --genome ${refGenome} --runMode ${runMode} ${"--blat " + blat_reject} ${"--pon " + m2_pon} ${"--gnomad " + gnomad} ${"--abra2 " + abra2_calls} 

    }

    runtime {
        docker         : "trinhanne/ngs_renv:0.1"       
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ~{diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
    File pass_judgement="${pairName}.judgement.vcf.gz"
    File pass_maf="${pairName}.judgement.maf"
    }
}
