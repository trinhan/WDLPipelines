## Author: Anne Trinh
## Last Modified: 23/11/2022
## This pipeline runs QCcheck, SNV, CNV and Manta
##
## Overall pipeline for the following
## - QC checks
## - SNV calling (consensus with Haplotypecaller, strelka2, pisces, vardict )
## - CNV calling
## - SV valling
##
## MODIFICATIONS:
##  1. Run Strelka2 somatic and germline
##  2. Merge VCF from all callers (Mutect 1/2, Strelka 2) prior to Oncotator and VEP. Annotated with INFO ME field: method used
##  3. Removed existing strelka implementation
##  4. Include option to run gatk4 and gatk3 
##  5. Test gatk4 outputs into absolute
##  6. include titan
##  7. consolidate gatk docker files and .jar files
##  8. allow hg38 compatibility
##  9. Include Manta and AnnotSV
##
## 
############################################
## DEFAULT OPTIONS (In the listed .json file)
#############################################
## refGenome = "hg38" (also works with "hg19")
## runQCCheck = false (Omit this step if picard metrics etc have already been run for this sample)
## targetedrun = false. Is a target panel used?

version 1.0

import "SNV-GermlineMulti.wdl" as SNV
import "Manta1.6.wdl" as Manta
import "annotsv.wdl" as AnnotSV
import "run_QC_checks.wdl" as runQC

workflow WGS_Germline_Workflow {
    input {
    File normalBam
    File normalBamIdx
    String pairName
    String ctrlName
    File refFasta
    File refFastaIdx
    File refFastaDict
    File targetIntervals
    File? pisces_reference
    String gatk_docker
    ## VEP input
    File vep_cache
    File? caddSnv
    File? caddSnvTbi
    File? caddIndel
    File? caddIndelTbi
    File? revel
    File? clinvar
    File? clinvarTbi
    String refGenome
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    String oncotree
    File AAlist
    String OncoKBtoken
    ## jointdiscovery inputs
    Array[File] HC_resources
    Array[File] HC_resources_index
    File gnomad
    File gnomadidx
    Int? minCallerSupport
    Int? HC_shard_counts
    ## Manta requirements
    File annotSVtar
    ## Other Options
    Boolean runQC
    Boolean targetedRun
    Boolean save_manta_evidence
    }

    String assembly = if refGenome=="hg19" then "GRCh37" else "GRCh38"
    String targName=basename(sub(targetIntervals,"\\.interval_list", ""))
    Int normalBam_size  = ceil(size(normalBam,  "G") + size(normalBamIdx,   "G")) 
    Int refFasta_size   = ceil(size(refFasta,   "G") + size(refFastaDict,   "G") + size(refFastaIdx, "G")) 
    Int db_snp_vcf_size = ceil(size(DB_SNP_VCF, "G")+size(DB_SNP_VCF_IDX, "G"))

    if (runQC){
        call runQC.PicardMultipleMetrics_Task as normalMM_Task {
            input:
                bam=normalBam,
                bamIndex=normalBamIdx,
                sampleName=ctrlName,
                refFasta=refFasta,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                targetIntervals=targetIntervals,
                baitIntervals=targetIntervals,
                gatk_docker=gatk_docker,
                refFasta_size=refFasta_size,
                db_snp_vcf_size=db_snp_vcf_size,
                bam_size=normalBam_size,
                targetedRun=targetedRun
        }
    }


    call SNV.runGermlineVariants as runSNV{
        input:
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            ctrlName=ctrlName,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            pisces_reference=pisces_reference,
            gatk_docker=gatk_docker,
            vep_cache=vep_cache,
            caddSnv=caddSnv,
            caddSnvTbi=caddSnvTbi,
            caddIndel=caddIndel,
            caddIndelTbi=caddIndelTbi,
            revel=revel,
            clinvar=clinvar,
            clinvarTbi=clinvarTbi,
            refGenome=refGenome,
            oncotree=oncotree,
            OncoKBtoken=OncoKBtoken,
            AAlist=AAlist,
            HC_resources=HC_resources,
            HC_resources_index=HC_resources_index,
            gnomad=gnomad,
            gnomadidx=gnomadidx,
            HC_shard_counts=HC_shard_counts,
            minCallerSupport=minCallerSupport
    }

    call Manta.MantaGermline as runManta {
        input: 
            sample_name = ctrlName,
            bam = normalBam,
            bam_index = normalBamIdx,
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIdx,
            disk_size=5*ceil(normalBam_size +refFasta_size),
            save_evidence=save_manta_evidence
    }

    call AnnotSV.annotsv as runAnnotSVManta {
        input:
          input_vcf=runManta.germline_sv_vcf ,
          input_vcf_idx=runManta.germline_sv_vcf_tbi ,
          genome_build=assembly,
          annotSVtar=annotSVtar,
          sampleName=ctrlName,
          typeofAnnotation="both",
          caller="Manta"
    }


output {
        ## Picard outputs
        File? Picard_QC_Output=normalMM_Task.picard_files
        File? Picard_HsMetrics=normalMM_Task.hsMetrics
        ## SNVoutputs
        File Merged_germline=runSNV.Merged_germline
        File Merged_germlineIdx=runSNV.Merged_germlineIdx
        File vep_annot = runSNV.vep_annot
        File? vep_summary_html=runSNV.vep_summary_html
        File oncokbMaf=runSNV.oncokbMaf
        ## MantaOutputs
        File manta_evidence_bam = runManta.evidence_bam
        File manta_evidence_bai = runManta.evidence_bai
        File mantaAnnotSV = runAnnotSVManta.sv_variants_tsv
    }

}