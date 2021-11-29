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

import "Pipeline_variantJudgment.wdl" as VCjudge
import "SNV-MultiCaller.wdl" as SomaticVC
import "abra2.wdl" as abra2
import "run_QC_checks.wdl" as runQC
import "cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl" as GATKCNVWorkflow
import "combine_tracks_version_modified.wdl" as CombineTracks
import "titan_workflow.wdl" as titan_workflow
import "VEP104.wdl" as VEP

workflow WGS_SNV_CNV_Workflow {
    input {
        # Bam files
        File tumorBam
        File? normalBam
        File tumorBamIdx
        File? normalBamIdx
        String pairName
        String? ctrlName
        String caseName
        # Ref Genome
        File refFasta
        File refFastaIdx
        File refFastaDict
        File targetIntervals
        ## Annotation Files
        File genome_bit
        File gnomad
        File gnomadidx
        File PoN
        File PoNidx    
        File? pisces_reference
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File DB_SNP_VCF
        File DB_SNP_VCF_IDX
        File? cosmicVCF
        File strelka_config
        ## QC Checkfiles
        File? regionFile
        File? captureNormalsDBRCLZip
        File? readGroupBlackList
        File HaplotypeDBForCrossCheck
        File? PicardMetrics_tumor
        ##CNV files required
        File centromere_tracks_seg
        File gistic_blacklist_tracks_seg
        File cnv_pon
        File common_snps
        ## VEP input
        File vep_cache
        File? caddSnv
        File? caddSnvTbi
        File? caddIndel
        File? caddIndelTbi
        File? revel
        File? clinvar
        File? clinvarTbi
        ##Other things required
        String gatk_docker
        String runMode
        String refGenome
        # Boolean optains
        Boolean run_CNQC
        Boolean forceComputePicardMetrics_tumor
        Boolean run_CrossCheck
        Boolean hasPicardMetrics_normal
        Boolean hasPicardMetrics_tumor
        Boolean forceComputePicardMetrics_normal
        Boolean run_functotar_cnv
        Boolean run_absolute
        Boolean run_VC_check
    }
    String assembly = if refGenome=="hg19" then "GRCh37" else "GRCh38"
    Int tumorBam_size=ceil(size(tumorBam, "G")+size(tumorBamIdx, "G"))

## Run the QC checks step here if specified
    call runQC.QCChecks as QCChecks {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            refFasta=refFasta,
            pairName=pairName,
            caseName=caseName,
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
            forceComputePicardMetrics_normal=forceComputePicardMetrics_normal,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            ctrlName=ctrlName
     }

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
            cosmicVCF=cosmicVCF,
            fracContam=QCChecks.fracContam,
            runMode=runMode,
            gatk_docker=gatk_docker,
            strelka_config=strelka_config
     }

    #### run the variant checking pipeline 
    if (run_VC_check) {
        call VCjudge.variantJudgement as VCjudge {
            input:
                merged_maf=somaticVC.Combined_raw_variants_maf,
                vcf=somaticVC.Combined_raw_variants_gz,
                vcfidx=somaticVC.Combined_raw_variants_tbi,
                indeltargets=somaticVC.VariantSitesBed,
                tumorBam=tumorBam,
                tumorBamIdx=tumorBamIdx,
                normalBam=normalBam,
                normalBamIdx=normalBamIdx,
                refFasta=refFasta,
                pairName=pairName,
                caseName=caseName,
                refFastaIdx=refFastaIdx,
                refFasta=refFasta,
                refFastaDict=refFastaDict,
                targetIntervals=targetIntervals,
                genome_bit=genome_bit,
                gnomad=gnomad,
                gnomadidx=gnomadidx,
                PoN=PoN,
                PoNidx=PoNidx,
                pisces_reference=pisces_reference,
                variants_for_contamination=variants_for_contamination,
                variants_for_contamination_idx=variants_for_contamination_idx,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                cosmicVCF=cosmicVCF, # can this be updated?
                fracContam=QCChecks.fracContam,
                runMode=runMode,
                gatk_docker=gatk_docker,
                strelka_config=strelka_config,
                refGenome=refGenome,
                tumorBam_size=tumorBam_size
        }
     }

    ## Change the output files depending on how it's done

    File new_variants_vcf = select_first([VCjudge.judgement_vcf, somaticVC.Combined_raw_variants_gz])
    
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

    if (defined(normalBam)==true && run_absolute) {
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
            gatk_docker=gatk_docker,
            group_id=pairName
        }
    # Estimate purity/ploidy, and from that compute absolute copy-number and mutation multiplicities.
        call absolute {
        input:
            maf=somaticVC.Combined_raw_variants_maf,
            seg_file=CombineTracksWorkflow.cnv_postprocessing_tumor_acs_seg,
            skew=CombineTracksWorkflow.cnv_postprocessing_tumor_acs_skew,
            pairName=pairName  
        }
    }    


    call titan_workflow.runtitan as runtitan {
    input:
        tumor_hets = GATK4CNV.het_allelic_counts_tumor,
        denoised_CR = GATK4CNV.denoised_copy_ratios_tumor,
        pairName = pairName,
        genomeBuild=refGenome
    }

    call VEP.variant_effect_predictor as vep {
        input:
          inputFile =new_variants_vcf,
          sample_name = pairName,
          assembly=assembly,
          species = "homo_sapiens",
          input_format = "vcf",
          cache=vep_cache,
          caddSnv=caddSnv,
          caddSnvTbi=caddSnvTbi,
          caddIndel=caddIndel,
          caddIndelTbi=caddIndelTbi,
          revelPlugin=revel,
          refFasta = refFasta,
          refFastaFai = refFastaIdx,
          clinvar=clinvar,
          clinvarTbi=clinvarTbi,
          gnomad=gnomad,
          gnomadIdx=gnomadidx
     }

    output {
        ####### QC check #######
        Float ContEst_contam = QCChecks.fracContam
        File? cross_check_fingprt_metrics=QCChecks.cross_check_fingprt_metrics
        File? copy_number_qc_report=QCChecks.copy_number_qc_report
        Array[File]? normal_picard_metrics=QCChecks.normal_bam_picard
        Array[File]? tumor_picard_metrics=QCChecks.tumor_bam_picard
        ####### SNV outputs ##########
        File? strelka2SomaticSNVs = somaticVC.strelka2SomaticSNVs
        File? strelka2SomaticIndels = somaticVC.strelka2SomaticIndels
        ####### pisces outputs ########
        File? pisces_tum_phased=somaticVC.pisces_tum_phased
        File? pisces_tum_unique=somaticVC.pisces_tum_unique
        File? pisces_venn=somaticVC.pisces_venn
        File? pisces_norm_same_site=somaticVC.pisces_norm_same_site
        File? pisces_tumor_variants=somaticVC.pisces_tumor_variants
        ####### M2 workflow2 outputs #####
        File M2_filtered_vcf=somaticVC.M2_filtered_vcf
        File M2_filtered_vcf_idx=somaticVC.M2_filtered_vcf_idx
        ####### merged output haplotypecaller
        File Combined_raw_variants_gz=somaticVC.Combined_raw_variants_gz
        File Combined_raw_variants_tbi=somaticVC.Combined_raw_variants_tbi
        File Combined_raw_variants_maf=somaticVC.Combined_raw_variants_maf
        File VariantSitesBed=somaticVC.VariantSitesBed
        ####### Copy Number - GATK4 CNV #######
        File? gatk4_het_allelic_counts_tumor = GATK4CNV.het_allelic_counts_tumor
        File? gatk4_normal_het_allelic_counts_tumor = GATK4CNV.het_allelic_counts_normal
        File? gatk4_copy_ratio_only_segments_tumor = GATK4CNV.copy_ratio_only_segments_tumor
        File? gatk4_modeled_segments_tumor = GATK4CNV.modeled_segments_tumor
        File? gatk4_denoised_copy_ratios_plot_tumor = GATK4CNV.denoised_copy_ratios_plot_tumor
        File? gatk4_modeled_segments_plot_tumor = GATK4CNV.modeled_segments_plot_tumor
        File? gatk4_called_copy_ratio_segments_tumor = GATK4CNV.called_copy_ratio_segments_tumor
        File gatk4_oncotated_called_file_tumor = select_first([GATK4CNV.oncotated_called_file_tumor, "null"])
        File gatk4_oncotated_called_gene_list_file_tumor = select_first([GATK4CNV.oncotated_called_gene_list_file_tumor, "null"])
        File gatk4_funcotated_called_file_tumor = select_first([GATK4CNV.funcotated_called_file_tumor, "null"])
        File gark4_funcotated_called_gene_list_file_tumor = select_first([GATK4CNV.funcotated_called_gene_list_file_tumor, "null"])
        ####### Absolute #######
        File? absolute_highres_plot=absolute.absolute_highres_plot
        File? absolute_rdata=absolute.absolute_rdata
        ###### Titan ##########
        File titan_optimal_params = runtitan.opt_params
        File titan_cluster_figs=runtitan.cluster_figures
        ###### VEP ##########
        File new_variants_gz=new_variants_vcf
        File vep_annot = vep.annotatedFile
        File? vep_summary_html=vep.summary_html
    }
}

task absolute {
    input {
    # TASK INPUT PARAMS
    File seg_file
    File maf
    String skew
    String pairName

    # RUNTIME INPUT PARAMS
    String preemptible ="3"
    String diskGB_boot ="14"
    String diskGB_buffer ="20"
    String machine_memoryGB = "7"
    String cpu ="1"
    }

    # COMPUTE DISK SIZE
    Int diskGB = ceil(size(seg_file, "G") + size(maf, "G") + diskGB_buffer)

    parameter_meta {
        seg_file : "filename pointing to the input - either a HAPSEG file or a segmentation file"
        maf : "filename pointing to a mutation annotation format (MAF) file. This specifies the data for somatic point mutations to be used by ABSOLUTE."
        skew : ""
        pairName : "a string for the name of the pair under analysis used for naming output files"
    }

    command {

        set -euxo pipefail

        SNV_MAF="${pairName}.snv.maf"
        INDEL_MAF="${pairName}.indel.maf"
        python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"
        
        grep -v "NA" ${seg_file} > no_nan_segs.tsv

        Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
        --seg_dat_fn no_nan_segs.tsv \
        --maf_fn $SNV_MAF \
        --indelmaf_fn $INDEL_MAF \
        --sample_name ${pairName} \
        --results_dir . \
        --ssnv_skew ${skew} \
        --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/
    }

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        # Plot showing the Purity/Ploidy values and the solutions
        File absolute_highres_plot="${pairName}.ABSOLUTE_plot.pdf"
        # An R file containing an object seg.dat which provides all of the information used to generate the plot.
        File absolute_rdata="${pairName}.PP-modes.data.RData"
    }
}

