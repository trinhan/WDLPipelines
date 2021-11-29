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
	Boolean runRealign
	# Annotation Files
	File genome_bit
	File gnomad
	File gnomadidx
	File PoN
	File PoNidx	
	File? pisces_reference
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
	}


	if (refGenome == "hg19"){
    call blat.blat as blat {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            MAF=merged_maf, 
            pairName=pairName,
            genome_bit=genome_bit,
            tumorBam_size=tumorBam_size,
            refGenome=refGenome
    }
	}

    if (runRealign) {
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
    		pisces_reference=pisces_reference,
    		variants_for_contamination=variants_for_contamination,
    		variants_for_contamination_idx=variants_for_contamination_idx,
    		DB_SNP_VCF=DB_SNP_VCF,
  			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
    		cosmicVCF=cosmicVCF, # can this be updated?
    		fracContam=fracContam,
		    runMode=runMode,
    		gatk_docker=gatk_docker,
    		strelka_config=strelka_config
    }
	}

	call PassJudgement {
		input:
			blat_reject=blat.rejectMaf,
			abra2_calls=somaticVC.Combined_raw_variants,
			vcf=vcf,
    		gnomad=gnomad,
    		m2_pon=PoN,
    		pairName=pairName
	}

	 output {
        File? blat_maf=blat.passMaf
        File? blat_reject=blat.rejectMaf
        File? abra2_tum_bam=runabra2.Abra2Tumour
        File? abra2_tum_bam_idx=runabra2.Abra2TumourIdx
        File? abra2_norm_bam=runabra2.Abra2Normal
        File? abra2_norm_bam_idx=runabra2.Abra2NormalIdx
        File? merged_abra2_calls=somaticVC.Combined_raw_variants
        File judgement_vcf=PassJudgement.pass_judgement
    }
}


task PassJudgement {
	input {
		File? blat_reject
		File? abra2_calls
		File vcf
    	File gnomad
    	File m2_pon
    	String pairName

    	# RUNTIME INPUT PARAMS
    	String? preemptible = "2"
    	String? diskGB_boot = "10"
    	String? diskGB_buffer ="5"
    	String? memoryGB ="4"
    	String? cpu ="1"
	}

	# DEFAULT VALUES
       
    Int diskGB = ceil(size(m2_pon, "G")+size(gnomad, "G"))*2  +  diskGB_buffer

	command {
		# Use an R script for this?
		# change blat to vcf file
		# bcftools isec for PoN
		
		Rscript /app/MergeVCFRecall.R --vcffile ${vcf} \
			--outputstr ${pairName} \
			--blat ${blat_reject} \ 
			--pon ${m2_pon} \
			${"--abra2 " + abra2_calls}
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
    }
}
