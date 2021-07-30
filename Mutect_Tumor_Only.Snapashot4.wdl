workflow CGA_Mutect_Tumor_Only {
	## WORKFLOW Description:
	## CGA workflow modified for no matched normal
	## Mutation only calls

    # WORKFLOW INPUT PARAMS
    # Configuration json file with optional parameters
    File cga_pipeline_config
    # Pair Input    
    # sample tumor BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File tumorBam
    # sample normal BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File tumorBamIdx
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    String caseName
    # a string for the name of the normal sample under analysis used for naming output files
    File readGroupBlackList
    # the FASTA file for the appropriate genome build (Reference sequence file)
    File refFasta
    # the FASTA file index for the reference genome (see http://www.htslib.org/doc/faidx.html)
    File refFastaIdx
    # the FASTA file dictionary for the reference genome (see https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)
    File refFastaDict
    # an interval list file that contains the locations of the targets
    File targetIntervals
    # an interval list file that contains the locations of the baits used
    File baitIntervals
    # VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs;
    # PROGRAMs whose CLP doesn't allow for this argument will quietly ignore it
    File DB_SNP_VCF
    # index file of VCF file of DB SNP variants
    File DB_SNP_VCF_IDX
    # catalogue of somatic mutations in VCF format    
    File cosmicVCF
    # 1000 genomes panel of normals in VCF format (MuTect1, MuTectFC)    
    File MuTectNormalPanel 
    # TSV file of chromsomal annotation ; chr, start, end, band, stain
    File cytoBandFile 
    # GenomeAnalysisToolkit jar file - a collection of command-line tools for analyzing high-throughput sequencing data
    File GATK4_JAR

    # Loading optional parameters from configuration file
    Map[String, String] runtime_params=read_json(cga_pipeline_config)

    # List of PONs for MAF filtering in task MafPonFilter
    File PONs_list
    Array[Object] PONs_data=read_objects(PONs_list)

    # COMPUTE FILE SIZE
    Int gatk4_jar_size  = ceil(size(GATK4_JAR,  "G"))
    Int tumorBam_size   = ceil(size(tumorBam,   "G") + size(tumorBamIdx,    "G")) 
    Int db_snp_vcf_size = ceil(size(DB_SNP_VCF, "G") + size(DB_SNP_VCF_IDX, "G"))
    Int refFasta_size   = ceil(size(refFasta,   "G") + size(refFastaDict,   "G") + size(refFastaIdx, "G"))

    # Does the sample already have picard metrics computed
    Boolean hasPicardMetrics_tumor           = false
    Boolean hasPicardMetrics_normal          = false
    # Should we compute picard metrics anyway, even if they exist
    Boolean forceComputePicardMetrics_tumor  = true
    Boolean forceComputePicardMetrics_normal = true
   
##############################


  call ContEST_Task {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            caseName=caseName,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            diskGB_buffer=runtime_params["ContEST_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["ContEST_Task.diskGB_boot"],
            preemptible=runtime_params["ContEST_Task.preemptible"],
            memoryGB=runtime_params["ContEST_Task.memoryGB"],
            cpu=runtime_params["ContEST_Task.cpu"]
    }


# Picard tasks (tumor and normal)
    ################################### 
    # The task runs 3 tools:
    # ValidateSamFile, CollectMultipleMetrics and CollectHSMetrics
    # ValidateSamFile makes sure the the given file is constructed correctly.
    # CollectMultipleMetrics collects multiple classes of metrics. This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time 
    # to cut down on the time spent reading in data from input files.
    # Available modules include CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, 
    # CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics.
    # CollectHSMetrics adds coverage statistics for WES files, on top of CollectMultipleMetrics.

    # tumor
    if (forceComputePicardMetrics_tumor || !hasPicardMetrics_tumor) {
        call PicardMultipleMetrics_Task as tumorMM_Task {
            input:
                bam=tumorBam,
                bamIndex=tumorBamIdx,
                sampleName=caseName,
                refFasta=refFasta,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                targetIntervals=targetIntervals,
                baitIntervals=baitIntervals,
                GATK4_JAR=GATK4_JAR,
                refFasta_size=refFasta_size,
                db_snp_vcf_size=db_snp_vcf_size,
                gatk4_jar_size=gatk4_jar_size,
                bam_size=tumorBam_size,
                validationStringencyLevel=runtime_params["PicardMultipleMetrics_Task.validationStringencyLevel"],
                run_clean_sam=runtime_params["PicardMultipleMetrics_Task.run_clean_sam"],
                diskGB_buffer=runtime_params["PicardMultipleMetrics_Task.diskGB_buffer"],
                diskGB_boot=runtime_params["PicardMultipleMetrics_Task.diskGB_boot"],
                preemptible=runtime_params["PicardMultipleMetrics_Task.preemptible"],
                memoryGB=runtime_params["PicardMultipleMetrics_Task.memoryGB"], 
                cpu=runtime_params["PicardMultipleMetrics_Task.cpu"]
        }
    }



    # PREPARE FOR SCATTER
    call CallSomaticMutations_Prepare_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals, # takes padded interval file (10bp on each side)
            nWay=runtime_params["CallSomaticMutations_Prepare_Task.nWay"],
            diskGB_boot=runtime_params["CallSomaticMutations_Prepare_Task.diskGB_boot"],
            preemptible=runtime_params["CallSomaticMutations_Prepare_Task.preemptible"]
    }


#SCATTER AND ANALYZE
    scatter (idx in CallSomaticMutations_Prepare_Task.scatterIndices) {
            # Identification of somatic point mutations in next generation sequencing data of cancer genomes.
            call Mutect1_Task {
                input:
                    tumorBam=tumorBam,
                    tumorBamIdx=tumorBamIdx,
                    caseName=caseName,
                    fracContam=ContEST_Task.fracContam,
                    mutectIntervals=CallSomaticMutations_Prepare_Task.interval_files[idx],
                    refFasta=refFasta,
                    refFastaIdx=refFastaIdx,
                    refFastaDict=refFastaDict,
                    DB_SNP_VCF=DB_SNP_VCF,
                    DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                    cosmicVCF=cosmicVCF,
                    readGroupBlackList=readGroupBlackList,
                    MuTectNormalPanel=MuTectNormalPanel,
                    refFasta_size=refFasta_size,
                    db_snp_vcf_size=db_snp_vcf_size,
                    tumorBam_size=tumorBam_size,
                    downsampleToCoverage=runtime_params["Mutect1_Task.downsampleToCoverage"],                    
                    diskGB_buffer=runtime_params["Mutect1_Task.diskGB_buffer"],
                    diskGB_boot=runtime_params["Mutect1_Task.diskGB_boot"],
                    preemptible=runtime_params["Mutect1_Task.preemptible"],
                    memoryGB=runtime_params["Mutect1_Task.memoryGB"],
                    cpu=runtime_params["Mutect1_Task.cpu"]
            }

            call Mutect2_Task {
                input:
                    tumorBam=tumorBam,
                    tumorBamIdx=tumorBamIdx,
                    caseName=caseName,
                    fracContam=ContEST_Task.fracContam,
                    mutectIntervals=CallSomaticMutations_Prepare_Task.interval_files[idx],
                    refFasta=refFasta,
                    refFastaIdx=refFastaIdx,
                    refFastaDict=refFastaDict,
                    readGroupBlackList=readGroupBlackList,
                    MuTectNormalPanel=MuTectNormalPanel,
                    GATK4_JAR=GATK4_JAR,
                    refFasta_size=refFasta_size,
                    tumorBam_size=tumorBam_size,
                    gatk4_jar_size=gatk4_jar_size,
                    diskGB_buffer=runtime_params["Mutect2_Task.diskGB_buffer"],
                    diskGB_boot=runtime_params["Mutect2_Task.diskGB_boot"], 
                    preemptible=runtime_params["Mutect2_Task.preemptible"],
                    memoryGB=runtime_params["Mutect2_Task.memoryGB"],
                    cpu=runtime_params["Mutect2_Task.cpu"]
            }
    }
			
 call MutectFC_Task {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            caseName=caseName,
            fracContam=ContEST_Task.fracContam,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            DB_SNP_VCF=DB_SNP_VCF,
            DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
            cosmicVCF=cosmicVCF,
            readGroupBlackList=readGroupBlackList,
            MuTectNormalPanel=MuTectNormalPanel,
            refFasta_size=refFasta_size,
            db_snp_vcf_size=db_snp_vcf_size,
            tumorBam_size=tumorBam_size,
            downsampleToCoverage=runtime_params["MutectFC_Task.downsampleToCoverage"],
            diskGB_buffer=runtime_params["MutectFC_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["MutectFC_Task.diskGB_boot"],
            preemptible=runtime_params["MutectFC_Task.preemptible"],
            memoryGB=runtime_params["MutectFC_Task.memoryGB"],
            cpu=runtime_params["MutectFC_Task.cpu"]
    }


  call GatherAndDeTiN_Task {
        input:
            caseName=caseName,
            mutect1_cs=Mutect1_Task.mutect1_cs,
            mutect2_cs=Mutect2_Task.mutect2_cs,
            release_version=runtime_params["GatherAndDeTiN_Task.release_version"],
            Mutation_prior=runtime_params["GatherAndDeTiN_Task.Mutation_prior"],
            TiN_prior=runtime_params["GatherAndDeTiN_Task.TiN_prior"],
            diskGB_buffer=runtime_params["GatherAndDeTiN_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["GatherAndDeTiN_Task.diskGB_boot"],
            preemptible=runtime_params["GatherAndDeTiN_Task.preemptible"],
            memoryGB=runtime_params["GatherAndDeTiN_Task.memoryGB"], 
            cpu=runtime_params["GatherAndDeTiN_Task.cpu"]
    }

# Oncotator is a tool for annotating human genomic point mutations and indels with data relevant to cancer researchers.
    call Oncotate_Task {
        input :
		    MUTECT1_CS=GatherAndDeTiN_Task.WXS_Mutation_M1_call_stats,
            MUTECT2_INDELS=GatherAndDeTiN_Task.WXS_Mutation_M2_call_stats,
            caseName=caseName,
            diskGB_buffer=runtime_params["Oncotate_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["Oncotate_Task.diskGB_boot"],
            preemptible=runtime_params["Oncotate_Task.preemptible"],
            memoryGB=runtime_params["Oncotate_Task.memoryGB"],
            cpu=runtime_params["Oncotate_Task.cpu"]
    }
    
  # VEP determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions. Simply input the coordinates of your variants and the nucleotide changes to find out the:
    call VEP_Task {
        input:
            mutect1_cs=GatherAndDeTiN_Task.WXS_Mutation_M1_call_stats,
            mutect2_vcf=GatherAndDeTiN_Task.WXS_Mutation_M2_call_stats,
            caseName=caseName,
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            refFasta_size=refFasta_size,
            diskGB_buffer=runtime_params["VEP_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["VEP_Task.diskGB_boot"],
            preemptible=runtime_params["VEP_Task.preemptible"],
            memoryGB=runtime_params["VEP_Task.memoryGB"],
            cpu=runtime_params["VEP_Task.cpu"]
    }


 call OrientationBias_filter_Task as oxoGOBF {
        input:   
            stub="oxog",
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            caseName=caseName,
            detailMetrics=tumorMM_Task.pre_adapter_detail_metrics,
            MAF=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
            GATK4_JAR=GATK4_JAR,
            refFasta=refFasta,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            gatk4_jar_size=gatk4_jar_size,
            diskGB_buffer=runtime_params["OrientationBias_filter_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["OrientationBias_filter_Task.diskGB_boot"],
            preemptible=runtime_params["OrientationBias_filter_Task.preemptible"],
            memoryGB=runtime_params["OrientationBias_filter_Task.memoryGB"],
            cpu=runtime_params["OrientationBias_filter_Task.cpu"]
    }

    # Detects and screens out FFPE artifacts from a set of SNV calls.
    # FFPE introduces multiple types of DNA damage including deamination, which converts cytosine to uracil and leads to downstream mispairing
    # in PCR: C>T/G>A. Because deamination occurs prior to ligation of palindromic Illumina adapters, likely deamination artifacts will have
    # a read orientation bias. The FFPE Filter Task uses this read orientation to identify artifacts and calculate a Phred scaled Q-score for FFPE artifacts.
    # .CG -> .TG <= DNA F1R2 (Context - ".CG", REF Allele - "C", ALT Allele - "T")
    # CG. -> CA.
    call OrientationBias_filter_Task as ffpeOBF {
        input:
            stub="ffpe",
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            caseName=caseName,
            detailMetrics=tumorMM_Task.pre_adapter_detail_metrics,
            MAF=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
            GATK4_JAR=GATK4_JAR,
            refFasta=refFasta,
            refFasta_size=refFasta_size,
            tumorBam_size=tumorBam_size,
            gatk4_jar_size=gatk4_jar_size,
            diskGB_buffer=runtime_params["OrientationBias_filter_Task.diskGB_buffer"],
            diskGB_boot=runtime_params["OrientationBias_filter_Task.diskGB_boot"],
            preemptible=runtime_params["OrientationBias_filter_Task.preemptible"],
            memoryGB=runtime_params["OrientationBias_filter_Task.memoryGB"],
            cpu=runtime_params["OrientationBias_filter_Task.cpu"]
    }

    # MAFPoNFilter uses a likelihood model to compare somatic mutations against a Panel of Normals (PoN)
    # in order to screen out somatic mutations. The PoN represents sequencing conditions in the case sample, 
    # including germline variants and technical artifacts. 

    scatter (pon_object in PONs_data) {
        call MAFPonFilter{
            input:
                MAFFile=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
                caseName=caseName,
                cytoBandFile=cytoBandFile,
                PONFile=pon_object.pon_url,
                stub=pon_object.pon_name,
                TOTNStr=runtime_params["MAFPonFilter.TOTNStr"],
                NMIN=runtime_params["MAFPonFilter.NMIN"],
                THRESH=runtime_params["MAFPonFilter.THRESH"],
                CODING_ONLY=runtime_params["MAFPonFilter.CODING_ONLY"],
                MIN_ALT_COUNT=runtime_params["MAFPonFilter.MIN_ALT_COUNT"],
                public_release=runtime_params["MAFPonFilter.public_release"],
                diskGB_buffer=runtime_params["MAFPonFilter.diskGB_buffer"],
                diskGB_boot=runtime_params["MAFPonFilter.diskGB_boot"],
                preemptible=runtime_params["MAFPonFilter.preemptible"],
                memoryGB=runtime_params["MAFPonFilter.memoryGB"],
                cpu=runtime_params["MAFPonFilter.cpu"]
        }
    }

    call blat {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            MAF=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf,
            caseName=caseName,
            tumorBam_size=tumorBam_size,
            diskGB_buffer=runtime_params["blat.diskGB_buffer"],
            diskGB_boot=runtime_params["blat.diskGB_boot"],
            preemptible=runtime_params["blat.preemptible"],
            memoryGB=runtime_params["blat.memoryGB"],
            cpu=runtime_params["blat.cpu"]
    }

    call merge_mafs_task {
        input:
            oxoGOBF_maf=oxoGOBF.allMaf,
            ffpeOBF_maf=ffpeOBF.allMaf,
            pon_filtered_mafs=MAFPonFilter.allMaf,
            blat_maf=blat.allMaf,
            caseName=caseName,
            diskGB_buffer=runtime_params["merge_mafs_task.diskGB_buffer"],
            diskGB_boot=runtime_params["merge_mafs_task.diskGB_boot"],
            preemptible=runtime_params["merge_mafs_task.preemptible"],
            memoryGB=runtime_params["merge_mafs_task.memoryGB"],
            cpu=runtime_params["merge_mafs_task.cpu"]
    }

 
  output {
        ####### Mutation Calling Tasks Outputs #######
        # MutectFC_Task
        File mutect_force_call_cs=MutectFC_Task.mutectfc_cs
        File mutect_force_call_pw=MutectFC_Task.mutectfc_pw
        File mutect_force_call_cw=MutectFC_Task.mutectfc_cw
        # Gathered MuTect1 and MuTect2 calls stats
        File mutect1_call_stats=GatherAndDeTiN_Task.WXS_Mutation_M1_call_stats
        File mutect2_call_stats=GatherAndDeTiN_Task.WXS_Mutation_M2_call_stats
        # deTiN outputs
        # Oncotator Output
        File mutect1_snv_mutect2_indel_maf=Oncotate_Task.WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf
        # Variant Effect
        File variant_effect_output=VEP_Task.VEP_Output
        File variant_effect_report=VEP_Task.VEP_Report
        # Merge MAF File Task
        File all_filters_passed_merged_maf=merge_mafs_task.merged_intersection_maf
        File after_filters_merged_maf=merge_mafs_task.merged_union_maf
        # Mutation Validator
        # File mutation_validator_pileup_preprocessing=mutation_validator.pileup_preprocessing_txt
        # File mutation_validator_validated_maf=mutation_validator.validated_maf
    
}}


# TASKS DEFINITION


task PicardMultipleMetrics_Task {

    # TASK INPUT PARAMS
    File bam
    File bamIndex
    String sampleName
    File targetIntervals
    File baitIntervals
    File refFasta
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    File GATK4_JAR

    String validationStringencyLevel
    String run_clean_sam

    # FILE SIZE
    Int bam_size
    Int refFasta_size
    Int gatk4_jar_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_stringencyLevel = "LENIENT" # may need to be adjusted depending on sequencing center
    String default_run_clean_sam = false

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(1.2*(ceil(bam_size + refFasta_size + gatk4_jar_size + db_snp_vcf_size + machine_diskGB_buffer)))

    String default_preemptible = if diskGB>100 then "0" else "1"   
    
    String stringencyLevel = if validationStringencyLevel != "" then validationStringencyLevel else default_stringencyLevel
    String clean_sam_flag = if run_clean_sam != "" then run_clean_sam else default_run_clean_sam

    parameter_meta {
        bam : "sample (normal or tumor) BAM file"
        bamIndex : "sample (normal or tumor) BAI file (BAM indexed)"
        sampleName : "sample (normal or tumor) name, prefix for output"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        DB_SNP_VCF_IDX : "dbSNP indexed file"
    }

    command <<<

        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} ValidateSamFile \
        --INPUT ${bam} \
        --OUTPUT "${sampleName}.bam_validation" \
        --MODE VERBOSE \
        --IGNORE_WARNINGS true \
        --REFERENCE_SEQUENCE ${refFasta} \
        --VALIDATION_STRINGENCY ${stringencyLevel}

        if [ "${clean_sam_flag}" = true ] ;
        then
            # Run bam through CleanSam to set MAPQ of unmapped reads to zero
            /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CleanSam \
            --INPUT ${bam} \
            --OUTPUT ${sampleName}.unmapped_reads_cleaned.bam
        fi

        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CollectMultipleMetrics \
        --INPUT ${bam} \
        --OUTPUT ${sampleName}.multiple_metrics \
        --REFERENCE_SEQUENCE ${refFasta} \
        --DB_SNP ${DB_SNP_VCF} \
        --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectInsertSizeMetrics \
        --PROGRAM QualityScoreDistribution \
        --PROGRAM MeanQualityByCycle \
        --PROGRAM CollectBaseDistributionByCycle \
        --PROGRAM CollectSequencingArtifactMetrics \
        --PROGRAM CollectQualityYieldMetrics \
        --PROGRAM CollectGcBiasMetrics
        
        #Extract OxoG metrics from generalized artifacts metrics. 
        # This tool extracts 8-oxoguanine (OxoG) artifact metrics from the output of CollectSequencingArtifactsMetrics 
        # (a tool that provides detailed information on a variety of
        # artifacts found in sequencing libraries) and converts them to the CollectOxoGMetrics tool's output format. This
        # conveniently eliminates the need to run CollectOxoGMetrics if we already ran CollectSequencingArtifactsMetrics in our
        # pipeline.
        /usr/local/jre1.8.0_73/bin/java -jar ${GATK4_JAR} ConvertSequencingArtifactToOxoG \
        --INPUT_BASE "${sampleName}.multiple_metrics" \
        --OUTPUT_BASE "${sampleName}.multiple_metrics.converted" \
        --REFERENCE_SEQUENCE ${refFasta} \
        --VALIDATION_STRINGENCY ${stringencyLevel}

        #zip up reports for QC Nozzle report
        zip ${sampleName}.picard_multiple_metrics.zip ${sampleName}.multiple_metrics.*

        # Collect WES HS metrics
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CollectHsMetrics \
        --INPUT ${bam} \
        --BAIT_INTERVALS ${targetIntervals} \
        --TARGET_INTERVALS ${baitIntervals} \
        --OUTPUT "${sampleName}.HSMetrics.txt" \
        --VALIDATION_STRINGENCY ${stringencyLevel}
    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File bam_validation="${sampleName}.bam_validation"
        File metricsReportsZip="${sampleName}.picard_multiple_metrics.zip"
        File alignment_summary_metrics="${sampleName}.multiple_metrics.alignment_summary_metrics"
        File bait_bias_detail_metrics="${sampleName}.multiple_metrics.bait_bias_detail_metrics"
        File bait_bias_summary_metrics="${sampleName}.multiple_metrics.bait_bias_summary_metrics"
        File base_distribution_by_cycle="${sampleName}.multiple_metrics.base_distribution_by_cycle.pdf"
        File base_distribution_by_cycle_metrics="${sampleName}.multiple_metrics.base_distribution_by_cycle_metrics"
        File gc_bias_detail_metrics="${sampleName}.multiple_metrics.gc_bias.detail_metrics"
        File gc_bias="${sampleName}.multiple_metrics.gc_bias.pdf"
        File gc_bias_summary_metrics="${sampleName}.multiple_metrics.gc_bias.summary_metrics"
        File insert_size_histogram="${sampleName}.multiple_metrics.insert_size_histogram.pdf"
        File insert_size_metrics="${sampleName}.multiple_metrics.insert_size_metrics"        
        File pre_adapter_detail_metrics="${sampleName}.multiple_metrics.pre_adapter_detail_metrics"
        File pre_adapter_summary_metrics="${sampleName}.multiple_metrics.pre_adapter_summary_metrics"
        File quality_by_cycle="${sampleName}.multiple_metrics.quality_by_cycle.pdf"
        File quality_by_cycle_metrics="${sampleName}.multiple_metrics.quality_by_cycle_metrics"
        File quality_distribution="${sampleName}.multiple_metrics.quality_distribution.pdf"
        File quality_distribution_metrics="${sampleName}.multiple_metrics.quality_distribution_metrics"
        File quality_yield_metrics="${sampleName}.multiple_metrics.quality_yield_metrics"
        File converted_oxog_metrics="${sampleName}.multiple_metrics.converted.oxog_metrics"
        File hsMetrics="${sampleName}.HSMetrics.txt"
    }
}


task ContEST_Task {
        
    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File refFasta
    File refFastaIdx
    File refFastaDict
    File targetIntervals    
    File SNP6Bed
    File HapMapVCF
    String caseName

    # FILE SIZE
    Int refFasta_size
    Int tumorBam_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer

    # COMPUTE DISK SIZE
    Int diskGB = ceil(tumorBam_size + refFasta_size
                + size(targetIntervals, "G") + size(SNP6Bed, "G") + size(HapMapVCF, "G")
                + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        caseName : "sample name, prefix for output"
        refFasta : "the FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        targetIntervals : ""
        SNP6Bed : ""
        HapMapVCF : "the population allele frequencies for each SNP in HapMap"
    }

    command <<<

        set -euxo pipefail
         
        java "-Xmx${command_memoryGB}g" -Djava.io.tmpdir=/tmp -jar /usr/local/bin/GenomeAnalysisTK.jar \
        -T ContEst \
        -I:eval,genotype ${tumorBam} \
        -L ${targetIntervals} \
        -L ${SNP6Bed} \
        -isr INTERSECTION \
        -R ${refFasta} \
        -l INFO \
        -pf ${HapMapVCF} \
        -o contamination.af.txt \
        --trim_fraction 0.03 \
        --beta_threshold 0.05 \
        -br contamination.base_report.txt \
        -mbc 100 \
        --min_genotype_depth 30 \
        --min_genotype_ratio 0.8

        python /usr/local/bin/extract_contamination.py contamination.af.txt fraction_contamination.txt \
        contamination_validation.array_free.txt ${caseName}

        #Contamination validation/consensus
        python /usr/local/populateConsensusContamination_v26/contaminationConsensus.py \
        --pass_snp_qc false \
        --output contest_validation.output.tsv \
        --annotation contamination_percentage_consensus_capture \
        --array contamination_validation.array_free.txt \
        --noarray contamination_validation.array_free.txt

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File contamDataFile="contamination_validation.array_free.txt"
        File contestAFFile="contamination.af.txt"
        File contestBaseReport="contamination.base_report.txt"
        File validationOutput="contest_validation.output.tsv"
        Float fracContam=read_float("fraction_contamination.txt")
    }
}


task CallSomaticMutations_Prepare_Task {

    # TASK INPUT PARAMS
    File targetIntervals
    File refFasta
    File refFastaIdx
    File refFastaDict

    String nWay

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot

    # DEFAULT VALUES
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_nWay = "10"

    String input_nWay = if nWay != "" then nWay else default_nWay

    parameter_meta {
        nWay : "Number of ways to scatter (MuTect1 and MuTect2)"
        mutectIntervals : "a list of genomic intervals over which MuTect1 will operate"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
    }

    command <<<

        set -euxo pipefail

        # Calculate disk size for all shards of the mutect run
        SIZE_FILE=split_base_sizes_disk.dat

        # Create list of indices for the scatter job
        seq 0 $((${input_nWay}-1)) > indices.dat

        # Run the prepare task that splits the .interval_list file into subfiles
        java -jar /usr/local/bin/GatkScatterGatherPrepare.jar . ${input_nWay} \
        --intervals ${targetIntervals} --reference_sequence ${refFasta}

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        memory         : "1 GB"
    }

    output {
        Array[File] interval_files=glob("gatk-scatter.*")
        Array[Int] scatterIndices=read_lines("indices.dat")
    }
}

task Mutect1_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    String caseName
    File mutectIntervals
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    File cosmicVCF
    File readGroupBlackList
    File MuTectNormalPanel
    File refFasta
    File refFastaIdx
    File refFastaDict
    Float fracContam

    String downsampleToCoverage

    # FILE SIZE
    Int tumorBam_size
    Int refFasta_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"    
    String default_downsampleToCoverage = "99999"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = 2*(ceil(tumorBam_size + refFasta_size + db_snp_vcf_size
                + size(mutectIntervals, "G") + size(cosmicVCF, "G") + size(readGroupBlackList, "G")
                + size(MuTectNormalPanel, "G") + select_first([diskGB_buffer, default_diskGB_buffer])))
    
    String downsample = if downsampleToCoverage != "" then downsampleToCoverage else default_downsampleToCoverage    

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        caseName : "a string for the name of the pair under analysis used for naming output files"
        mutectIntervals : "a list of genomic intervals over which MuTect1 will operate"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        cosmicVCF : "catalogue of somatic mutations in VCF format"
        readGroupBlackList : ""
        MuTectNormalPanel : "1000 genomes panel of normals in VCF format"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "fraction of cross-sample contamination, output from ContEst task"       
        downsampleToCoverage : "downsample reads to a given capping threshold coverage"
    }

    command <<<

        set -euxo pipefail

        #variable for normal panel
        NORMAL_PANEL_FLAG_AND_VAL=""
        if [ -s "${MuTectNormalPanel}" ] ; then
            NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${MuTectNormalPanel}" ;
        fi ;

        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
        -L ${mutectIntervals} \
        --tumor_sample_name ${caseName} \
        -I:tumor ${tumorBam} \
        --reference_sequence ${refFasta} \
        --fraction_contamination ${fracContam} \
        --dbsnp ${DB_SNP_VCF} \
        --cosmic ${cosmicVCF} \
        --read_group_black_list ${readGroupBlackList} \
        --out ${caseName}.MuTect1.call_stats.txt \
        --coverage_file ${caseName}.MuTect1.coverage.wig.txt \
        --power_file ${caseName}.MuTect1.power.wig.txt \
        --downsample_to_coverage ${downsample} \
        $NORMAL_PANEL_FLAG_AND_VAL

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File mutect1_cs="${caseName}.MuTect1.call_stats.txt"
        File mutect1_pw="${caseName}.MuTect1.power.wig.txt"
        File mutect1_cw="${caseName}.MuTect1.coverage.wig.txt"
    }
}


task Mutect2_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    String caseName
    File mutectIntervals
    File readGroupBlackList
    File MuTectNormalPanel
    File refFasta
    File refFastaIdx
    File refFastaDict
    Float fracContam
    File GATK4_JAR

    # FILE SIZE
    Int tumorBam_size
    Int refFasta_size
    Int gatk4_jar_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = 2*(ceil(tumorBam_size +  refFasta_size + gatk4_jar_size
                + size(mutectIntervals, "G") + size(readGroupBlackList, "G") + size(MuTectNormalPanel, "G")
                + machine_diskGB_buffer))

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        caseName : "tumor sample name, prefix for output"
        mutectIntervals : "a list of genomic intervals over which MuTect2 will operate" 
        DB_SNP_VCF : ""
        readGroupBlackList : ""
        MuTectNormalPanel : "1000 genomes panel of normals in VCF format"
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "fraction of cross-sample contamination, output from ContEst task"
        GATK4_JAR : ""
    }

    command <<<

        set -euxo pipefail

        #variable for normal panel
        NORMAL_PANEL_FLAG_AND_VAL=""
        if [ -s "${MuTectNormalPanel}" ] ; then
            BZ="${MuTectNormalPanel}.gz"
            #bgzip the file and index it
            bgzip ${MuTectNormalPanel} # 0 level of compression flag -- compress without compression
            tabix $BZ 
            NORMAL_PANEL_FLAG_AND_VAL="--normal_panel $BZ" ;
        fi ;

        #MuTect2 wants names that match those in the BAMs so grab them from the BAMs
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} GetSampleName \
        -I ${tumorBam} -O tumorName.txt
         TUMOR_NAME=`cat tumorName.txt`

        #mutect 2 ----- gatk4
        /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} Mutect2 \
        --input ${tumorBam} \
        --tumor-sample "$TUMOR_NAME" \
        --reference ${refFasta} \
        --panel-of-normals $BZ \
        --contamination-fraction-to-filter ${fracContam} \
        --intervals ${mutectIntervals} \
        --output "${caseName}.MuTect2.call_stats.unfiltered.unaf.txt"

        #filter the variants
        /usr/local/jre1.8.0_73/bin/java -jar -Xmx4g ${GATK4_JAR} FilterMutectCalls \
        -O ${caseName}.MuTect2.call_stats.filtered.unaf.txt -V ${caseName}.MuTect2.call_stats.unfiltered.unaf.txt
		# not sure what this does
       # python /usr/local/bin/process_af.py "${caseName}.MuTect2.call_stats.filtered.unaf.txt"  \
       # "${caseName}.MuTect2.call_stats.txt" "$TUMOR_NAME" "Normal" "${caseName}" "Normal"
		
    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File mutect2_cs="${caseName}.MuTect2.call_stats.filtered.unaf.txt"
   #     File process_af_file="/usr/local/bin/process_af.py"
     }
}


task MutectFC_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    String caseName
    File mutectIntervals
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    File cosmicVCF
    File readGroupBlackList
    File MuTectNormalPanel
    File refFasta
    File refFastaIdx
    File refFastaDict
    Float fracContam

    String downsampleToCoverage

    # FILE SIZE
    Int tumorBam_size
     Int refFasta_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_downsampleToCoverage = "99999"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + refFasta_size + db_snp_vcf_size
                + size(mutectIntervals, "G") + size(cosmicVCF, "G") + size(readGroupBlackList, "G")
                + size(MuTectNormalPanel, "G") + machine_diskGB_buffer)

    String downsample = if downsampleToCoverage != "" then downsampleToCoverage else default_downsampleToCoverage

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
         caseName : "tumor sample name, prefix for output"
        mutectIntervals : ""
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        cosmicVCF : ""
        readGroupBlackList : ""
        MuTectNormalPanel : "1000 genomes panel of normals in VCF format"
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "fraction of cross-sample contamination, output from ContEst task"
        downsampleToCoverage : ""
    }

    command <<<
        
        set -euxo pipefail

        #variable for normal panel
        NORMAL_PANEL_FLAG_AND_VAL=""
        if [ -s "${MuTectNormalPanel}" ] ; then
            NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${MuTectNormalPanel}" ;
        fi ;

        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
        -L ${mutectIntervals} \
        --tumor_sample_name ${caseName} \
        -I:tumor ${tumorBam}  \
        --reference_sequence ${refFasta} \
        --fraction_contamination ${fracContam} \
        --dbsnp ${DB_SNP_VCF} \
        --cosmic ${cosmicVCF} \
        --force_output \
        --read_group_black_list ${readGroupBlackList} \
        --out ${caseName}.MuTect1.call_stats.txt \
        --coverage_file ${caseName}.MuTect1.coverage.wig.txt \
        --power_file ${caseName}.MuTect1.power.wig.txt \
        --downsample_to_coverage ${downsample} \
        $NORMAL_PANEL_FLAG_AND_VAL

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File mutectfc_cs="${caseName}.MuTect1.call_stats.txt"
        File mutectfc_pw="${caseName}.MuTect1.power.wig.txt"
        File mutectfc_cw="${caseName}.MuTect1.coverage.wig.txt"
    }
}

task GatherAndDeTiN_Task {

    # TASK INPUT PARAMS
    Array[File] mutect1_cs
    Array[File] mutect2_cs
    String caseName
    File? seg_file
    File? tumor_hets
    File? normal_hets
    File? exac_pickle

    String TiN_prior
    String Mutation_prior
    String release_version

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_TiN_prior = "0.2"
    String default_Mutation_prior = "0.05"
    String default_release_version = "v1.8.5"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(seg_file, "G") + size(tumor_hets, "G") + size(normal_hets, "G")
                    + size(exac_pickle, "G") + machine_diskGB_buffer)

    String Prior_TiN = if TiN_prior != "" then TiN_prior else default_TiN_prior
    String Prior_Mutation = if Mutation_prior != "" then Mutation_prior else default_Mutation_prior
    String deTiN_version = if release_version != "" then release_version else default_release_version

    parameter_meta {
        mutect1_cs : ""
        mutect2_cs : ""
        caseName: "a string for the name of the pair under analysis used for naming output files"
        seg_file : "filename pointing to the input - either a HAPSEG file or a segmentation file"
        tumor_hets : "a tab delimited file with read counts from the tumor sample."
        normal_hets : "a tab delimited file with read counts from the normal sample."
        exac_pickle : ""
        TiN_prior : ""
        Mutation_prior : ""
    }

    command <<<

        set -euxo pipefail

python <<CODE

mutect1_cs_file_array = '${sep="," mutect1_cs}'.split(",")
mutect2_cs_file_array = '${sep="," mutect2_cs}'.split(",")
mutect1_cs_list = open('mutect1_cs_list.txt', 'w')
mutect2_cs_list = open('mutect2_cs_list.txt', 'w')

for i in range(len(mutect1_cs_file_array)):
    mutect1_cs_list.write(mutect1_cs_file_array[i] + '\n')
for i in range(len(mutect2_cs_file_array)):
    mutect2_cs_list.write(mutect2_cs_file_array[i] + '\n')
mutect1_cs_list.close()
mutect2_cs_list.close()
CODE

        ######## Files #################

        # MuTect1 Files
        MUTECT1_CS="${caseName}.MuTect1.call_stats.txt"
        # MuTect2 Files
        MUTECT2_CS="${caseName}.MuTect2.call_stats.vcf"
        MUTECT2_INDELS="${caseName}.MuTect2.call_stats.indels.vcf"
        MUTECT2_OTHER="${caseName}.MuTect2.call_stats.other.vcf"

        ######## Merge MuTect1 and MuTect2 Call Stats #################

        # Gather MuTect1 call stats
        python /usr/local/bin/merge_callstats.py "mutect1_cs_list.txt" $MUTECT1_CS
        # Gather MuTect2 call stats
        python /usr/local/bin/merge_callstats.py "mutect2_cs_list.txt" $MUTECT2_CS

        # Filter out indels and SNPs for MuTect2 results
        python /usr/local/bin/filter_indels.py $MUTECT2_CS $MUTECT2_INDELS $MUTECT2_OTHER

        ######## Pull deTiN from GitHub ###########################        
    #     git clone https://github.com/broadinstitute/deTiN.git
#         git -C deTiN/ checkout tags/${deTiN_version}

        ######## Run deTiN on MuTect 1 & 2 ###########################
#         generate null outputs
#         echo null > ${caseName}.TiN_estimate.txt
#         echo null > ${caseName}.TiN_estimate_CI.txt
#         echo null > ${caseName}.deTiN_SSNVs.txt
#         echo null > ${caseName}.deTiN_indels.txt
#         echo null > ${caseName}.deTiN_aSCNAs.txt
#         echo null > ${caseName}_SSNVs_plot.png
#         echo null > ${caseName}_KmeansEval_plot.png
#         echo null > ${caseName}_TiN_models_plot.png
#         echo null > ${caseName}_TiN_hets_aSCNA_model.png
#         echo null > ${caseName}_KmeansEval_scatter_plot.png
# 
#         python deTiN/deTiN/deTiN.py --mutation_data_path $MUTECT1_CS --cn_data_path ${seg_file} \
#         --tumor_het_data ${tumor_hets} --normal_het_data ${normal_hets} --exac_data_path ${exac_pickle} \
#         --output_name ${caseName} --TiN_prior ${Prior_TiN} --mutation_prior ${Prior_Mutation} \
#         --indel_data_path $MUTECT2_INDELS --indel_data_type "MuTect2"

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File WXS_Mutation_M1_call_stats="${caseName}.MuTect1.call_stats.txt"
        File WXS_Mutation_M2_call_stats="${caseName}.MuTect2.call_stats.vcf"
        File WXS_Mutation_M2_indels_only="${caseName}.MuTect2.call_stats.indels.vcf"
      #   Float TiN=read_float("${caseName}.TiN_estimate.txt")
#         Int number_added_SSNV=read_int("${caseName}.number_of_SSNVs_added.txt")
#         String TiN_CI=read_string("${caseName}.TiN_estimate_CI.txt")
#         File deTiN_call_stats="${caseName}.deTiN_SSNVs.txt"
#         File deTiN_SSNVs_plot="${caseName}_SSNVs_plot.png"
#         File deTiN_aSCNA_kmeans_RSS_plot="${caseName}_KmeansEval_plot.png"
#         File deTiN_aSCNA_scatter_plot="${caseName}_KmeansEval_scatter_plot.png"
#         File deTiN_TiN_modes_plot="${caseName}_TiN_models_plot.png"
#         File deTiN_indels="${caseName}.deTiN_indels.txt"
#         File aSCNA_model="${caseName}_TiN_hets_aSCNA_model.png"
#         File deTiN_segments="${caseName}.deTiN_aSCNAs.txt"
    }
}

task Oncotate_Task {

    # TASK INPUT PARAMS
    File MUTECT1_CS
    File MUTECT2_INDELS
    String caseName
    File oncoDBTarBall

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "15"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
   
    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(oncoDBTarBall, "G") * 5 + machine_diskGB_buffer)

    parameter_meta {
        MUTECT1_CS : ""
        MUTECT2_INDELS : ""
        caseName : "tumor sample name, prefix for output"
          oncoDBTarBall : ""
    }

    command <<<

        set -x

        ######## Files #################

        # MuTect1 files
        MUTECT1_CS_PASSED="${caseName}.MuTect1.call_stats.passed.txt"
        MUTECT1_CS_REJECTED="${caseName}.MuTect1.call_stats.rejected.txt"
        MUTECT1_CS_MAFLITE="${caseName}.MuTect1.call_stats.maflite"
        MUTECT1_CS_ANNOTATED_MAF="${caseName}.MuTect1.call_stats.maflite.annotated.maf"

        # MuTect2 files
        MUTECT2_INDELS_PASSED="${caseName}.MuTect2.call_stats.passed.vcf"
        MUTECT2_INDELS_REJECTED="${caseName}.MuTect2.call_stats.rejected.vcf"
        MUTECT2_INDELS_ANNOTATED_MAF="${caseName}.MuTect2.call_stats.indels.annotated.maf"
        
      
        # Merged files
        MUTECT2_STRELKA_MERGED_MAF="${caseName}.MuTect2_INDEL.Strelka_INDEL.merged.maf" 
        MUTECT1_MUTECT2_STRELKA_MERGED_MAF="${caseName}.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.maf"
        MUTECT1_MUTECT2_STRELKA_ANNOTATED_VCF="${caseName}.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.vcf"

        ######## Processing call stats before annotation #################

        # Filter MuTect1 mutation calls that passed filter
        python /usr/local/bin/filter_passed_mutations.py ${MUTECT1_CS} $MUTECT1_CS_PASSED $MUTECT1_CS_REJECTED "KEEP"
        # Convert MuTect1 call stats to MafLite
        python /usr/local/bin/callstats_to_maflite.py $MUTECT1_CS_PASSED $MUTECT1_CS_MAFLITE

        # Filter MuTect2 mutation calls that passed filter
        python /usr/local/bin/filter_passed_mutations.py ${MUTECT2_INDELS} $MUTECT2_INDELS_PASSED $MUTECT2_INDELS_REJECTED "PASS"

             ########################### Unzip Oncotator database ###############################

        #find TARBALL type 
        TYPE=`echo 'if("${oncoDBTarBall}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ;

        #obtain the name of the directory for oncodb and unpack based on the TYPE
        if [ "$TYPE" == "GZ" ] ; then
            ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ;
            tar -xzf ${oncoDBTarBall}
        else
            ONCO_DB_DIR_NAME=`tar -tf ${oncoDBTarBall} |head -1` ;
            tar -xf ${oncoDBTarBall} ;
        fi ;

        ################## Annotate MuTect1, MuTect2, Strelka call stats ###########################

        # Annotate MuTect1 call stats (MAFLITE to TCGAMAF)
        # --infer-onps \
        /root/oncotator_venv/bin/oncotator -i MAFLITE --db-dir `pwd`/$ONCO_DB_DIR_NAME \
        --longer-other-tx \
        -a tumor_barcode:${caseName} \
        $MUTECT1_CS_MAFLITE $MUTECT1_CS_ANNOTATED_MAF hg19

        # Annotate MuTect2 call stats (VCF to TCGAMAF)
        /root/oncotator_venv/bin/oncotator -i VCF --db-dir `pwd`/$ONCO_DB_DIR_NAME \
        --longer-other-tx \
        --skip-no-alt \
        -a tumor_barcode:${caseName} \
        $MUTECT2_INDELS_PASSED $MUTECT2_INDELS_ANNOTATED_MAF hg19

      
        ####################### Merge MuTect1, MuTect2, Strelka annotated TCGAMAFs into one file ####################

        # Merge Strelka and MuTect2 indels (only Strelka calls selected, MuTect2 and Strelka common calls are annotated)
      #  python /usr/local/bin/merge_strelka_mutect2.py $STRELKA_ANNOTATED_MAF $MUTECT2_INDELS_ANNOTATED_MAF $MUTECT2_STRELKA_MERGED_MAF
        # Merge MuTect1 SNVs and Strelka & MuTect2 indels into one MAF file
        python /usr/local/bin/maf_merge.py $MUTECT1_CS_ANNOTATED_MAF $MUTECT2_INDELS_ANNOTATED_MAF $MUTECT1_MUTECT2_STRELKA_MERGED_MAF

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File WXS_Mutation_M1_SNV_M2_INDEL_Strelka_INDEL_annotated_maf="${caseName}.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.maf"               
    }
}


task VEP_Task {

    # TASK INPUT PARAMS
    File mutect1_cs
    File mutect2_vcf
    String caseName
    File VEP_File
    File refFasta
    File refFastaIdx
    File refFastaDict

    # FILE SIZE
    Int refFasta_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(mutect2_vcf, "G") + size(mutect1_cs, "G") + size(VEP_File, "G") * 5 + refFasta_size + machine_diskGB_buffer)

    parameter_meta {
        mutect2_vcf : ""
        mutect1_cs : ""
         caseName : "tumor sample name, prefix for output"
        VEP_File : ""
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
    }

    command <<<

        set -euxo pipefail

        ######################## Files #####################################

        MUTECT1_VCF="${caseName}.MuTect1.call_stats.vcf"
        MERGED_VCF="${caseName}.MuTect1.MuTect2.Strelka.vcf"

        ############## Pre-process MuTect1 call stats #################################

        # convert MuTect1 call stats into VCF
        python /usr/local/bin/callstats_to_vcf.py ${mutect1_cs} $MUTECT1_VCF "${caseName}" "Normal"

        ################ Merge MuTect1, MuTect2 and Strelka VCFs into one ##############################
        # combine variants into one VCF file
        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/GenomeAnalysisTK.jar -T CombineVariants \
        -R ${refFasta} \
        --variant $MUTECT1_VCF \
        --variant ${mutect2_vcf} \
        --out $MERGED_VCF \
        --assumeIdenticalSamples  \
        -U ALLOW_SEQ_DICT_INCOMPATIBILITY

        ############### Run Variant Effector Predictor (VEP)  #########################

        #make a link from the home directory to the current directory to avoid running out of disk space on the boot disk
        mkdir -v vep_data_dir
        #delete the existing directory first to make a successful link
        rm -rf ~/.vep
        ln -vs `pwd`/vep_data_dir ~/.vep

        #Run the merged RAW VCF from the call stats through VEP
        #In either case unpack the data into the home directory where VEP expects to find it
        IS_ZIP=`echo ${VEP_File}|grep -Pic '\.zip$'` ;
        if [ "$IS_ZIP" -eq "1" ] ;
        then
            #it's a zip file
            unzip -d ~ ${VEP_File} ;
        else
            #tar ball
            tar -C ~ -xvzf ${VEP_File}
        fi ;
        /ensembl-tools-release-83/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
        --vcf -i $MERGED_VCF --symbol --cache --vcf --offline 

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {        
        File VEP_Output="variant_effect_output.txt"
        File VEP_Report="variant_effect_output.txt_summary.html"
    }
}


task OrientationBias_filter_Task {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File? detailMetrics
    File MAF
    String caseName
    String stub
    File refFasta
    File GATK4_JAR

    # FILE SIZE
    Int tumorBam_size
    Int refFasta_size
    Int gatk4_jar_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "7"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + refFasta_size + gatk4_jar_size + size(MAF, "G") + machine_diskGB_buffer)

    parameter_meta {
        tumorBam: "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM)"
        MAF : "filename pointing to a mutation annotation format (MAF) file (data for somatic point mutations)"
        detailMetrics : "output from the Picard Multiple Metrics CollectSequencingArtifactMetrics run with the settings here; this allows for passthrough instead of recomputing"       
        refFasta : "Reference that was used to align BAM"
        stub : "string used to indicate in the output the effect name"
    }

    command <<<

        set -euxo pipefail

        SNV_MAF="${caseName}.snv.maf"
        INDEL_MAF="${caseName}.indel.maf"
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"

        ################################

        python /usr/local/bin/get_context_ref_alt_alleles.py -s ${stub} -i ${caseName}

        CONTEXT=$( cat "${caseName}.context.${stub}.txt" )
        ARTIFACT_ALLELE=$( cat "${caseName}.artifact_allele.${stub}.txt" )

        if [ -s "${detailMetrics}" ] ;
        then
            DETAIL_METRICS_FILE=${detailMetrics}
        else
            /usr/local/jre1.8.0_73/bin/java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CollectSequencingArtifactMetrics \
            --INPUT ${tumorBam} \
            --OUTPUT ${caseName} \
            --REFERENCE_SEQUENCE ${refFasta}
            DETAIL_METRICS_FILE="${caseName}.pre_adapter_detail_metrics"  
        fi ;

        # Now parsing metrics for Q value
        python /usr/local/orientationBiasFilter/annotate_orientationBiasQ.py -i ${caseName} -m $DETAIL_METRICS_FILE -c $CONTEXT -a $ARTIFACT_ALLELE

        Q=$( cat "${caseName}.orientation_BiasQ.txt" )

        # Appending Q value to MAF
        python /usr/local/orientationBiasFilter/AppendAnnotation2MAF.py -i ${caseName} -m $SNV_MAF -f ${stub}_Q -v $Q

        OUTPUT_INTERVAL_FILE="${caseName}.intervals"
        OUTPUT_INFO_FILE="${caseName}.orientation_info.txt"

        python /usr/local/orientationBiasFilter/write_interval_file.py "${caseName}.${stub}_Q.maf.annotated" $OUTPUT_INTERVAL_FILE

        java "-Xmx${command_memoryGB}g" -jar /usr/local/orientationBiasFilter/GenomeAnalysisTK.jar \
        --analysis_type OxoGMetrics \
        -R ${refFasta} \
        -I ${tumorBam} \
        -L $OUTPUT_INTERVAL_FILE \
        -o $OUTPUT_INFO_FILE

        # Now appending orientation bias information to the MAF
        python /usr/local/orientationBiasFilter/AppendOrientationBiasFields2MAF.py \
        -i ${caseName} -m "${caseName}.${stub}_Q.maf.annotated" -b ${tumorBam} \
        -f $OUTPUT_INFO_FILE

        REF_ALLELE_COMP=$( cat "${caseName}.ref_allele_compliment.${stub}.txt" )
        ARTIFACT_ALLELE_COMP=$( cat "${caseName}.artifact_allele_compliment.${stub}.txt")

        bash -c "source /matlab_source_file_2012a.sh && /usr/local/orientationBiasFilter/orientationBiasFilter ${caseName}.OrientationBiasInfo.maf \
        ${caseName}.OrientationBiasFilter.maf . '0' '1' '0.96' '0.01' '-1' '30' '1.5' \
        $REF_ALLELE_COMP $ARTIFACT_ALLELE_COMP i_${stub}" ;

        #########################################

        #merge back indels into OBF output
        python /usr/local/bin/tsvConcatFiles.py $INDEL_MAF "${caseName}.OrientationBiasFilter.maf" \
        --outputFilename="${caseName}.OrientationBiasFilter.${stub}.indel_snp_merged.filtered.maf"

        python /usr/local/bin/tsvConcatFiles.py $INDEL_MAF "${caseName}.OrientationBiasFilter.unfiltered.maf" \
        --outputFilename="${caseName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.maf"

        python /usr/local/bin/add_judgement_column.py \
        --input "${caseName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.maf" \
        --output "${caseName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.with_judgement_annotations.maf" \
        --column "i_${stub}_cut" \
        --pass_flag "0"
        
        zip -r ${caseName}.${stub}_OBF_figures.zip ./figures
        zip ${caseName}.${stub}.maf.zip *.maf

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File zipped_mafs="${caseName}.${stub}.maf.zip"
        Float q_val=read_float("${caseName}.orientation_BiasQ.txt")
        File OBF_figures="${caseName}.${stub}_OBF_figures.zip"
        Int num_passed_mutations=read_int("${caseName}.OrientationBiasFilter.maf.pass_count.txt")
        Int num_rejected_mutations=read_int("${caseName}.OrientationBiasFilter.maf.reject_count.txt")
        File passMaf="${caseName}.OrientationBiasFilter.${stub}.indel_snp_merged.filtered.maf"        
        File allMaf="${caseName}.OrientationBiasFilter.${stub}.indel_snp_merged.unfiltered.with_judgement_annotations.maf"
    }
}

task MAFPonFilter {

    # TASK INPUT PARAMS
    File MAFFile
    String caseName
    File PONFile
    File cytoBandFile
    String stub

    String TOTNStr
    String NMIN
    String THRESH
    String CODING_ONLY
    String MIN_ALT_COUNT
    String public_release

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_NMIN = "1"
    String default_THRESH = "-2.5"
    String default_CODING_ONLY = "0"
    String default_MIN_ALT_COUNT = "3"
    String default_TOTNStr = "TN"
    String default_public_release = "false"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(MAFFile, "G") + size(PONFile, "G") + size(cytoBandFile, "G") + machine_diskGB_buffer)

    String input_TOTNStr = if TOTNStr != "" then TOTNStr else default_TOTNStr
    String input_NMIN = if NMIN != "" then NMIN else default_NMIN
    String input_THRESH = if THRESH != "" then THRESH else default_THRESH
    String input_CODING_ONLY = if CODING_ONLY != "" then CODING_ONLY else default_CODING_ONLY
    String input_MIN_ALT_COUNT  = if MIN_ALT_COUNT != "" then MIN_ALT_COUNT else default_MIN_ALT_COUNT
    String input_public_release = if public_release != "" then public_release else default_public_release

    parameter_meta {
        MAFFile : "input MAF for PON filter analysis"
        PONFile : "formatted panel-of-normals file"
        cytoBandFile : "TSV file of chromosomal annotation: chr, start, end, band, stain"
        caseName : "used to name the output files and title other output"
        TOTNStr : "indicating pair status : can be 'TO' or 'TN'"
        NMIN : "minimum count of samples in a token PoN bin that are used to estimate the PoN log-likelihood (to reduce the effect of tumor contamination in the PoN)"
        thresh : "threshold for pass/fail with log likelihood"
        WCUT : "threshold for pass/fail with pon-computed weight"
        CODING_ONLY : "analyze coding regions only?"
        MIN_ALT_COUNT : "the minimum t_alt_count for mutations in the maf that pass the maf_pon_filter"
    }

    command <<<

        set -euxo pipefail
 
        #the cytoBand file should be in the expected location
        mkdir -pv /xchip/cga/reference/annotation/db/ucsc/hg19/
        cp -v  ${cytoBandFile} /xchip/cga/reference/annotation/db/ucsc/hg19/cytoBand.txt
        ls -alht /xchip/cga/reference/annotation/db/ucsc/hg19/cytoBand.txt

        #run the filter
        #Note "." is used for the parameter file whose usage is *not* exposed/functional in this WDL
        bash -c "source /matlab_source_file_2013a.sh && /usr/local/bin/maf_pon_filter \
        ${MAFFile} ${PONFile} ${caseName} ${input_TOTNStr} . ${input_NMIN} ${input_THRESH} \
        0.5 . ${input_CODING_ONLY} ${input_MIN_ALT_COUNT}" ;

        # Count number of passed and rejected mutations
        python /usr/local/bin/count_variants.py "${caseName}.pon_annotated.pass.maf" "${caseName}.count_passed_mutations.txt"
        python /usr/local/bin/count_variants.py "${caseName}.pon.blacklist.txt" "${caseName}.count_rejected_mutations.txt"

        python /usr/local/bin/add_judgement_column.py \
        --input "${caseName}.pon_annotated.maf" \
        --output "${caseName}.${stub}.with_judgement.all.maf" \
        --column "pon_pass_loglike" \
        --pass_flag "1"

        if [ "${input_public_release}" = true ] ; 
        then
            python /usr/local/bin/remove_columns.py \
            --IN_MAF "${caseName}.pon_annotated.pass.maf" \
            --OUT_MAF "${caseName}.${stub}.pass.maf"

            python /usr/local/bin/remove_columns.py \
            --IN_MAF "${caseName}.${stub}.with_judgement.all.maf" \
            --OUT_MAF "${caseName}.${stub}.all.maf"

        else

            python /usr/local/bin/add_stub.py \
            --input "${caseName}.pon_annotated.pass.maf" \
            --output "${caseName}.${stub}.pass.maf" \
            --stub "${stub}"

            python /usr/local/bin/add_stub.py \
            --input "${caseName}.${stub}.with_judgement.all.maf" \
            --output "${caseName}.${stub}.all.maf" \
            --stub "${stub}"

        fi

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb: if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        Int num_rejected_mutations=read_int("${caseName}.count_rejected_mutations.txt")
        Int num_passed_mutations=read_int("${caseName}.count_passed_mutations.txt")
        File allMaf="${caseName}.${stub}.all.maf"
        File passMaf="${caseName}.${stub}.pass.maf"
    }
}


task blat {

    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File MAF
    File hg19_bit
    String caseName    

    # FILE SIZE
    Int tumorBam_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + size(MAF, "G") + size(hg19_bit, "G") + machine_diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM)"
        MAF : "filename pointing to a mutation annotation format (MAF) file (data for somatic point mutations)"
        caseName : "tumor sample name, string prefix of the output"
    }

    command <<<

        set -euxo pipefail

        cp -v ${hg19_bit} /opt/hg19.2bit
        python /opt/realign.py ${tumorBam} ${MAF} ${caseName}

        # Count number of passed and rejected mutations
        python /usr/local/bin/count_variants.py "${caseName}.blat.maf" "${caseName}.count_passed_mutations.txt"
        python /usr/local/bin/count_variants.py "${caseName}.blat.rejected.maf" "${caseName}.count_rejected_mutations.txt"

        python /usr/local/bin/add_judgement_column.py \
        --input "${caseName}.blat.all.maf" \
        --output "${caseName}.blat.all.with_judgement_annotations.maf" \
        --column "realign_judgment" \
        --pass_flag "KEEP"       

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        Int num_rejected_mutations=read_int("${caseName}.count_rejected_mutations.txt")
        Int num_passed_mutations=read_int("${caseName}.count_passed_mutations.txt")
        File passMaf="${caseName}.blat.maf"
        File debug_results="${caseName}.blat.rejected.maf"        
        File allMaf="${caseName}.blat.all.with_judgement_annotations.maf"
    }
}

task merge_mafs_task {

    # TASK INPUT PARAMS
    File oxoGOBF_maf
    File ffpeOBF_maf
    Array[File] pon_filtered_mafs   
    File blat_maf
    String caseName

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(size(oxoGOBF_maf, "G") + size(ffpeOBF_maf, "G") + size(blat_maf, "G") + machine_diskGB_buffer)

    parameter_meta {
        oxoGOBF_maf       : "MAF with annotations of passing or failing OxoG Orientation Bias filter"
        ffpeOBF_maf       : "MAF with annotations of passing or failing FFPE Orientation Bias filter"        
        blat_maf          : "MAF with annotations of passing or failing BLAT Re-aligner filter"
        pon_filtered_mafs : "MAF(s) with annotations of passing or failing MAFPONFilter"
    }

    command <<<

        set -euxo pipefail

python <<CODE

pon_filtered_mafs_array = '${sep="," pon_filtered_mafs}'.split(",")

command_string_file = open('command_string.txt', 'w')
command_string = ""
for filename in pon_filtered_mafs_array: 
    pon_name = filename.split(".")[1]
    command_string = command_string + "--{0}={1} ".format(pon_name, filename)

command_string_file.write(command_string)
command_string_file.close()
CODE

        COMMAND_STRING=`cat command_string.txt`

        python /usr/local/bin/merge_mafs_after_filters.py \
        --pair_name=${caseName} \
        --oxogOBF=${oxoGOBF_maf} \
        --ffpeOBF=${ffpeOBF_maf} \
        --blat=${blat_maf} \
        $COMMAND_STRING

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File merged_intersection_maf="${caseName}.passed_all_filters.maf"
        File merged_union_maf="${caseName}.merged_union.maf"
    }
}


task mutation_validator {

    # TASK INPUT PARAMS
    String caseName
    File MAF
    File tumorBam
    File tumorBamIdx
    File? tumorRNABam
    File? tumorRNABamIdx

    String maf_type

    # FILE SIZE
    Int tumorBam_size
  
    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "10"
    String default_preemptible = "1"
    String default_diskGB_boot = "15"
    String default_diskGB_buffer = "20"
    String default_maf_type = "wex"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB

    # COMPUTE DISK SIZE
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer
    Int diskGB = ceil(tumorBam_size + size(MAF, "G") 
                        + size(tumorRNABam, "G") + size(tumorRNABamIdx, "G")
                        + machine_diskGB_buffer)

    String input_maf_type = if maf_type != "" then maf_type else default_maf_type

    parameter_meta {
        tumorBam     : "sample tumor BAM file"
        tumorBamIdx  : "sample tumor BAI file (indexed BAM file)"
        caseName     : "a string for the name of the pair under analysis used for naming output files"

    }

    command <<<

        set -euxo pipefail

        SNV_MAF="${caseName}.snv.maf"
        INDEL_MAF="${caseName}.indel.maf"
        PREPROCESSED_FILE="${caseName}.pileup_preprocessing.txt"
        VALIDATED_SNV_MAF="${caseName}.snp.validated.maf"
        VALIDATED_INDEL_MAF="${caseName}.indel.validated.maf"
        VALIDATED_MAF="${caseName}.validated.maf"

        #RNATYPE="hg19-chr"
        RNATYPE="hg19"

        # split MAF into INDELs and SNPs
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${MAF} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"

        python /opt/src/fh_MutationValidatorPreprocess/validation_wrapper_firehose_library_hack.py \
        --mafsnp $SNV_MAF \
        --mafindel $INDEL_MAF \
        --wextumor ${tumorBam} \
        --rnatumor ${default="None" tumorRNABam} \
        --rnatype $RNATYPE \
        --out ${caseName} 

        bash -c "source /matlab_source_file_2012a.sh && /opt/src/fh_MutationValidator/mutation_validator_wrapper $PREPROCESSED_FILE $SNV_MAF ${input_maf_type} ${caseName}.snp"
        bash -c "source /matlab_source_file_2012a.sh && /opt/src/fh_MutationValidator/mutation_validator_wrapper $PREPROCESSED_FILE $INDEL_MAF ${input_maf_type} ${caseName}.indel"
        
        python /usr/local/bin/tsvConcatFiles.py $VALIDATED_SNV_MAF $VALIDATED_INDEL_MAF \
        --outputFilename=$VALIDATED_MAF

    >>>

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : if diskGB_boot != "" then diskGB_boot else default_diskGB_boot
        preemptible    : if preemptible != "" then preemptible else default_preemptible
        cpu            : if cpu != "" then cpu else default_cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }

    output {
        File pileup_preprocessing_txt="${caseName}.pileup_preprocessing.txt"
        File validated_maf="${caseName}.validated.maf"
    }
}