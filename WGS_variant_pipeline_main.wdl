## Author: Anne Trinh
## Last Modified: 18/11/2022
## This pipeline is a modification of the Getz somatic variant pipeline modified for hg38
## Legacy features have been retained in this pipeline and can be turned off using Boolean input
##
## Overall pipeline for the following
## - QC checks
## - SNV calling
## - CNV calling
## - SV valling
## - judgement (Indel realigner)
## - TITAN and/or ABSOLUTE
## - vep of the final outputs
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
## runMode = "Paired" (can also be "TumOnly")
## runQCCheck = false (Omit this step if picard metrics etc have already been run for this sample)
## run_CNQC = false. Run copy number QC in QC step?
## run_CrossCheck = false. Run Cross Check Fingerprints?
## run_CNQC = false. Run copy number QC in QC step?
## run_CrossCheck = false. Run Cross Check Fingerprints?
## run_blat = false. Run blat (hg19 only)?
## run_abra2 = false. Run the indel realigner (abra2)?
## run_absolute = false. Run Absolute to determine ploidy/purity?
## run_functotar_cnv = false. Annotae the CNVs with funcotator?
## run_SV_paired = false. Run paired SV-calling (no point if matched normal is blood)
## targetedrun = false. Is a target panel used?

version 1.0

import "Pipeline_variantJudgment.wdl" as VCjudge
import "SNV-MultiCaller.wdl" as SomaticVC
import "abra2.wdl" as abra2
import "run_QC_checks.wdl" as runQC
import "cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl" as GATKCNVWorkflow
import "combine_tracks_version_modified.wdl" as CombineTracks
import "titan_workflow.wdl" as titan_workflow
import "VEP104.wdl" as VEP
import "Manta1.6.wdl" as Manta
import "annotsv.wdl" as AnnotSV

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
        File? genome_bit
        File gnomad
        File gnomadidx
        File PoN
        File PoNidx
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
        File cnv_intervals
        ## SV region files
        File SVRegionBed
        File SVRegionBedTbi
        File annotSVtar
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
        Boolean runQCCheck
        Boolean run_CNQC
        Boolean run_Picard_tumor
        Boolean run_CrossCheck
        Boolean run_Picard_normal
        Boolean run_functotar_cnv
        Boolean run_absolute
        Boolean run_blat
        Boolean run_abra2
        Boolean targetedRun
        Boolean run_SV_paired
        Boolean save_manta_evidence=true
        Float? fracContam =0.01 # by default set at 0.01
    }
    String assembly = if refGenome=="hg19" then "GRCh37" else "GRCh38"
    Int tumorBam_size=ceil(size(tumorBam, "G")+size(tumorBamIdx, "G"))
    Int normalBam_size=ceil(size(normalBam, "G")+size(normalBamIdx, "G"))
    Boolean run_VC_check = if (run_blat || run_abra2) then true else false

## Run the QC checks step here if specified
    if (runQCCheck){
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
            run_Picard_tumor=run_Picard_tumor, 
            run_CNQC=run_CNQC,
            run_CrossCheck=run_CrossCheck,
            run_Picard_normal=run_Picard_normal,
            normalBam=normalBam,
            fracContam=fracContam,
            normalBamIdx=normalBamIdx,
            ctrlName=ctrlName,
            targetedRun=targetedRun
     }
    }

    Float new_contam_frac = select_first([fracContam,QCChecks.fracContam])

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
            variants_for_contamination=variants_for_contamination,
            variants_for_contamination_idx=variants_for_contamination_idx,
            DB_SNP_VCF=DB_SNP_VCF,
            DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
            cosmicVCF=cosmicVCF,
            fracContam=new_contam_frac,
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
                ctrlName=ctrlName,
                refFastaIdx=refFastaIdx,
                refFasta=refFasta,
                refFastaDict=refFastaDict,
                targetIntervals=targetIntervals,
                genome_bit=genome_bit,
                gnomad=gnomad,
                gnomadidx=gnomadidx,
                PoN=PoN,
                PoNidx=PoNidx,
                variants_for_contamination=variants_for_contamination,
                variants_for_contamination_idx=variants_for_contamination_idx,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                cosmicVCF=cosmicVCF, # can this be updated?
                fracContam=new_contam_frac,
                runMode=runMode,
                gatk_docker=gatk_docker,
                strelka_config=strelka_config,
                refGenome=refGenome,
                tumorBam_size=tumorBam_size,
                runAbra2Realign=run_abra2,
                runBlatRealign=run_blat
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
            intervals=cnv_intervals,
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
          gnomadIdx=gnomadidx,
          canonical=true,
          hgvs=true,
          protein=true,
          polyphen="b",
          af_1kg=true,
          af_gnomad=true,
          ccds=true,
          domains=true,
          pubmed=true,
          regulatory=true,
          symbol=true,
          uniprot=true,
          biotype=true
     }

    call Manta.MantaSomaticSV as MantaWF {
        input:
          sample_name = caseName,
          tumor_bam = tumorBam,
          tumor_bam_index = tumorBamIdx,
          normal_bam = if (run_SV_paired) then normalBam else "null",
          normal_bam_index = if (run_SV_paired) then normalBamIdx else "null",
          ref_fasta = refFasta,
          ref_fasta_index = refFastaIdx,
          region_bed=SVRegionBed,
          disk_size=5*ceil(tumorBam_size + (if (run_SV_paired) then normalBam_size else 0)),
          region_bed_index=SVRegionBedTbi,
          save_evidence=save_manta_evidence
    } 

    call AnnotSV.annotsv as AnnotSVWF {
        input:
          input_vcf=MantaWF.somatic_sv_vcf,
          input_vcf_idx=MantaWF.somatic_sv_vcf_tbi,
          genome_build=assembly,
          annotSVtar=annotSVtar,
          sampleName=caseName,
          typeofAnnotation="both",
          caller="Manta"
    }


    output {
        ####### QC check #######
        Float? ContEst_contam = QCChecks.fracContam
        File? cross_check_fingprt_metrics=QCChecks.cross_check_fingprt_metrics
        File? copy_number_qc_report=QCChecks.copy_number_qc_report
        File? normal_picard_metrics=QCChecks.normal_bam_picard
        File? tumor_picard_metrics=QCChecks.tumor_bam_picard
        File? tumor_hsmetrics=QCChecks.tumor_bam_hsmetrics
        File? normal_hsmetrics=QCChecks.normal_bam_hsmetrics
        File? tumor_cleaned_unmapped_bam=QCChecks.tumor_cleaned_unmapped_bam
        File? normal_cleaned_unmapped_bam=QCChecks.normal_cleaned_unmapped_bam
        ####### SNV outputs ##########
        File? strelka2SomaticSNVs = somaticVC.strelka2SomaticSNVs
        File? strelka2SomaticIndels = somaticVC.strelka2SomaticIndels
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
        ##### Manta outputs #####
        File manta = AnnotSVWF.sv_variants_tsv
        File? manta_evidence_bam = MantaWF.evidence_bam
        File? manta_evidence_bai = MantaWF.evidence_bai
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
    Int diskGB = ceil(size(seg_file, "G")) + ceil(size(maf, "G")) + select_first([diskGB_buffer, 10])

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
        disks          : "local-disk " + diskGB + " HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        # Plot showing the Purity/Ploidy values and the solutions
        File absolute_highres_plot="${pairName}.ABSOLUTE_plot.pdf"
        # An R file containing an object seg.dat which provides all of the information used to generate the plot.
        File absolute_rdata="${pairName}.PP-modes.data.RData"
    }
}

