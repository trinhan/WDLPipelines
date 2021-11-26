## Overall pipeline for the following
## - QC checks
## - SNV calling
## - CNV calling
## - judgement
## TITAN and/or ABSOLUTE
## vep and oncokb of the final outputs
##
## MODIFICATIONS:
##  1. Allow optional run for absolute and cnv
##  2. Run Strelka2 somatic and germline
##  3. Merge VCF from all callers (Mutect 1/2, Strelka 2) prior to Oncotator and VEP. Annotated with INFO ME field: method used
##  4. Removed existing strelka implementation
##  5. Include option to run gatk4 and gatk3 
##  6. Test gatk4 outputs into absolute
##  7. Run funcotator instead of oncotator
##  8. include titan
##  9. consolidate gatk docker files and .jar files
##  10. allow hg38 compatibility
##
## TODO MODIFICATIONS:
##  4. include varscan

version 1.0

import "blat_runner.wdl" as blat
import "SNVMultiCaller" as SomaticVC
import "abra2.wdl" as abra2
import "run_QC_checks.wdl" as runQC
import "https://github.com/broadinstitute/gatk/blob/4.1.7.0/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl" as GATKCNVWorkflow
import "combine_tracks_version_modified.wdl" as CombineTracks

workflow WGS_SNV_CNV_Workflow {
	input {
		# Bam files
		File tumorBam
 	    File normalBam
 	    File? tumorBamIdx
	    File? normalBamIdx
	    String pairName
		String? ctrlName
		String caseName
		# Ref Genome
		File refFasta
		File refFastaIdx
		File refFastaDict
		File targetIntervals
    	String ctrlName
    	## Annotation Files
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
		## Annotation for CNV

		# CNV files required
    	File centromere_tracks_seg
    	File gistic_blacklist_tracks_seg
    	# Other things required
    	String gatk_docker
    	String runMode
    	# Boolean optains

	}


## Run the QC checks step here if specified
	call runQC.QCChecks as QCChecks{
		input:
    		normalBam=normalBam,
        	normalBamIdx=normalBamIdx,
        	tumorBam=tumorBam,
        	tumorBamIdx=tumorBamIdx,
        	refFasta=refFasta,
    		pairName=pairName,
 			caseName=caseName,
 			ctrlName=ctrlName,
    		refFastaIdx=refFastaIdx,
    		refFastaDict=refFastaDict,
			targetIntervals=targetIntervals,
			baitIntervals=targetIntervals,
    		gatk_docker=gatk_docker,
    		refGenome=refGenome,
    		DB_SNP_VCF=DB_SNP_VCF,
  			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
    		gnomad=gnomad,
    		gnomad_idx=gnomadidx,
    		captureNormalsDBRCLZip=captureNormalsDBRCLZip,
    		regionFile=regionFile,
    		readGroupBlackList=readGroupBlackList,
    		HaplotypeDBForCrossCheck=HaplotypeDBForCrossCheck, 
    		run_CNQC=run_CNQC,
    		run_CrossCheck=run_CrossCheck,
    		hasPicardMetrics_tumor=hasPicardMetrics_tumor,
    		hasPicardMetrics_normal=hasPicardMetrics_normal,
    		forceComputePicardMetrics_tumor=forceComputePicardMetrics_tumor,
    		forceComputePicardMetrics_normal=forceComputePicardMetrics_normal
    }

 File picardMetrics=select_first([ tumorMM_Task.pre_adapter_detail_metrics, PicardMetrics_tumor])

# Call somatic variant calling
    call SomaticVC.runVariantCallers as somaticVC {
    	input:
    		normalBam=normalBam,
        	normalBamIdx=normalBamIdx,
        	tumorBam=tumorBam,
        	tumorBamIdx=tumorBamIdx,
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
    		fracContam= QCChecks.fracContam,
		    runMode=runMode,
    		gatk_docker=gatk_docker,
    		strelka_config=strelka_config
    }
	
# Call CNV 
       call GATKCNVWorkflow.CNVSomaticPairWorkflow as GATK4CNV {
        input:
            tumor_bam=tumorBam,
            tumor_bam_idx=tumorBamIdx,
            normal_bam=normalBam,
            normal_bam_idx=normalBamIdx,
            ref_fasta=refFasta,
            ref_fasta_fai=refFastaIdx,
            ref_fasta_dict=refFastaDict,
            intervals=targetIntervals,
            common_sites=common_snps,
            read_count_pon=cnv_pon,
            is_run_funcotator=run_functotar_cnv,
            gatk_docker=gatk_docker
        }

    if (defined(normalBam)==true){
        call CombineTracks.CombineTracksWorkflow as CombineTracksWorkflow {
        input:
            tumor_called_seg = GATK4CNV.called_copy_ratio_segments_tumor,
            tumor_modeled_seg = GATK4CNV.modeled_segments_tumor,
            af_param = GATK4CNV.allele_fraction_parameters_tumor,
            matched_normal_called_seg = select_first([GATK4CNV.called_copy_ratio_segments_normal, "null"]),
            ref_fasta=refFasta,
            ref_fasta_fai=refFastaIdx,
            ref_fasta_dict=refFastaDict,
            centromere_tracks_seg=centromere_tracks_seg,
            gistic_blacklist_tracks_seg =gistic_blacklist_tracks_seg,
            columns_of_interest = columns_of_interest,
            group_id=pairName
        }
    }

}