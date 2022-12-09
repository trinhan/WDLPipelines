## Author: Anne Trinh
## Last Modified: 23/11/2022
## This pipeline runs QCcheck, SNV, CNV and Manta
##
## Overall pipeline for the following
## - QC checks
## - SNV calling (consensus with Haplotypecaller, strelka2, pisces, vardict )
## - CNV calling (GATK, using an existing panel of normals)
## - SV valling
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
import "cnv_wdl/germline/cnv_germline_case_workflow.wdl" as CNV

workflow WGS_Germline_Workflow {
    input {
    File normalBam
    File normalBamIdx
    String ctrlName
    File refFasta
    File refFastaIdx
    File refFastaDict
    File snvtargetIntervals
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
    ## CNV requirements
    File cnv_contig_ploidy_model_tar
    File cnv_filtered_intervals
    Array[File]+ gcnv_model_tars
    File cnvIntervals
    Int cnv_max_number_events
    Int cnv_max_number_events_passed
    Int cnv_num_intervals_per_Scatter
    Int ref_copy_number_autosomal_contigs
    File cnv_blacklist_intervals
    Int gcnv_max_copy_number

    ## Other Options
    Boolean runQC
    Boolean targetedRun
    Boolean save_manta_evidence
    }

    String assembly = if refGenome=="hg19" then "GRCh37" else "GRCh38"
    String targName=basename(sub(snvtargetIntervals,"\\.interval_list", ""))
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
                targetIntervals=snvtargetIntervals,
                baitIntervals=snvtargetIntervals,
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
            targetIntervals=snvtargetIntervals,
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

    call CNV.CNVGermlineCaseWorkflow as runCNV {
        input:
            intervals=cnvIntervals,
            blacklist_intervals=cnv_blacklist_intervals,
            filtered_intervals=cnv_filtered_intervals,
            normal_bams=normalBam,
            normal_bais=normalBamIdx,
            sample_name=ctrlName,
            contig_ploidy_model_tar=cnv_contig_ploidy_model_tar,
            gcnv_model_tars=gcnv_model_tars,
            num_intervals_per_scatter=50000,
            ref_fasta_dict=refFastaDict,
            ref_fasta_fai=refFastaIdx,
            ref_fasta=refFasta,
            gatk_docker=gatk_docker,
            ref_copy_number_autosomal_contigs=ref_copy_number_autosomal_contigs,
            allosomal_contigs=["chrX","chrY"],
            maximum_number_events=cnv_max_number_events,
            maximum_number_pass_events=cnv_max_number_events_passed,
            gcnv_max_copy_number=gcnv_max_copy_number
    }

    call AnnotSV.annotsv as runAnnotSVCNV {
        input:
          input_vcf=runCNV.genotyped_segments_vcfs,
          input_vcf_idx=runCNV.genotyped_segments_vcf_indexes,
          genome_build=assembly,
          annotSVtar=annotSVtar,
          sampleName=ctrlName,
          typeofAnnotation="both",
          caller="gCNV"
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
        File? manta_evidence_bam = runManta.evidence_bam
        File? manta_evidence_bai = runManta.evidence_bai
        File mantaAnnotSV = runAnnotSVManta.sv_variants_tsv
        ## CNV outputs
        File denoised_copy_ratios = runCNV.denoised_copy_ratios
        File genotyped_intervals_vcf_indexes =runCNV.genotyped_intervals_vcf_indexes
        File genotyped_intervals_vcfs = runCNV.genotyped_intervals_vcfs
        File genotyped_segments_vcf_indexes=runCNV.genotyped_segments_vcf_indexes
        File genotyped_segments_vcfs=runCNV.genotyped_segments_vcfs
        File sample_contig_ploidy_calls_tars=runCNV.sample_contig_ploidy_calls_tars
        File gCNV_annotSV = runAnnotSVCNV.sv_variants_tsv
    }

}