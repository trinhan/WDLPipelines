version 1.0
## quick mod: override cross check lane for hg19 only
workflow QCChecks {
    input {
    File tumorBam
    # sample normal BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File? normalBam
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File tumorBamIdx
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File? normalBamIdx
    # a string for the name of the pair under analysis used for naming output files
    String pairName
    # a string for the name of the tumor sample under analysis used for naming output files
    String caseName
    # a string for the name of the normal sample under analysis used for naming output files
    String? ctrlName
    # list of read groups to exclude from the analysis in MuTect1 and MuTect_FC tasks
    File refFasta
    File refFastaIdx
    File refFastaDict
    File targetIntervals
    File baitIntervals
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX
    String gatk_docker
    File gnomad
    File gnomad_idx
    # For CNQC
    File? captureNormalsDBRCLZip
    File? regionFile
    File? readGroupBlackList

    ## Option: refGenome "hg19" or "hg18"
    String refGenome
    # ContEst
    # CrossCheckLane
    File HaplotypeDBForCrossCheck
    Float? fracContam

    # Does the sample already have picard metrics computed
    Boolean hasPicardMetrics_tumor = false
    Boolean hasPicardMetrics_normal = false
    # Forceto compute picard metrics anyway, even if they exist
    Boolean forceComputePicardMetrics_tumor = true
    Boolean run_CNQC = false
    Boolean run_CrossCheck = true
    Boolean forceComputePicardMetrics_normal = if defined (normalBam) then true else false
    # Is the run WGS (false) or a targeted hybrid capture panel (true)
    Boolean targetedRun
    }

    Int tumorBam_size   = ceil(size(tumorBam,   "G") + size(tumorBamIdx,    "G")) 
    Int normalBam_size  = if defined(normalBam) then ceil(size(normalBam,  "G") + size(normalBamIdx,   "G")) else 0
    Int db_snp_vcf_size = ceil(size(DB_SNP_VCF, "G") + size(DB_SNP_VCF_IDX, "G"))
    Int refFasta_size   = ceil(size(refFasta,   "G") + size(refFastaDict,   "G") + size(refFastaIdx, "G"))

    File capNorm=select_first([captureNormalsDBRCLZip, "null"])
    File regF=select_first([regionFile, "null"])
    File readGrpBL=select_first([readGroupBlackList, "null"])

    Boolean override_CNQC = if defined (normalBam) then run_CNQC else false
    Boolean runContEst = if defined(fracContam) || !defined(normalBam) then false else true


    if (override_CNQC){
        call CopyNumberReportQC_Task {
            input:
                tumorBam=tumorBam,
                tumorBamIdx=tumorBamIdx,
                normalBam=normalBam,
                normalBamIdx=normalBamIdx,
                readGroupBlackList=readGrpBL,
                tumorBam_size=tumorBam_size,
                normalBam_size=normalBam_size,
                captureNormalsDBRCLZip=capNorm,
                regionFile=regF
        }
    }
    if (runContEst){
	call ContEST_Task {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            targetIntervals=targetIntervals,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            gatk_docker=gatk_docker
    }
}

    # Program to check that all read groups within the set of BAM files appear to come from the same individual.
    if (run_CrossCheck){
        call CrossCheckLaneFingerprints_Task {
            input:
                tumorBam=tumorBam,
                normalBam=normalBam,
                tumorBamIdx=tumorBamIdx,
                normalBamIdx=normalBamIdx,
                pairName=pairName,
                tumorBam_size=tumorBam_size,
                normalBam_size=normalBam_size,
                HaplotypeDBForCrossCheck=HaplotypeDBForCrossCheck,
                gatk_docker=gatk_docker
        }
    }

    #Picard tasks (tumor and normal)
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
                gatk_docker=gatk_docker,
                refFasta_size=refFasta_size,
                db_snp_vcf_size=db_snp_vcf_size,
                bam_size=tumorBam_size,
                targetedRun=targetedRun
        }
    }

    #normal
    if (forceComputePicardMetrics_normal || !hasPicardMetrics_normal && defined(normalBam)) {
        call PicardMultipleMetrics_Task as normalMM_Task {
            input:
                bam=normalBam,
                bamIndex=normalBamIdx,
                sampleName=ctrlName,
                refFasta=refFasta,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                targetIntervals=targetIntervals,
                baitIntervals=baitIntervals,
                gatk_docker=gatk_docker,
                refFasta_size=refFasta_size,
                db_snp_vcf_size=db_snp_vcf_size,
                bam_size=normalBam_size,
                targetedRun=targetedRun
        }
    }

    output {

    	####### QC Tasks Outputs #######
        # Copy Number QC Report files
        File? tumor_bam_lane_list=CopyNumberReportQC_Task.tumorBamLaneList
        File? normal_bam_lane_list=CopyNumberReportQC_Task.normalBamLaneList
        File? tumor_bam_read_coverage_lane=CopyNumberReportQC_Task.tumorRCL
        File? normal_bam_read_coverage_lane=CopyNumberReportQC_Task.normalRCL
        File? copy_number_qc_report=CopyNumberReportQC_Task.CopyNumQCReport
        File? copy_number_qc_report_png=CopyNumberReportQC_Task.CopyNumQCReportPNG
        Int? copy_number_qc_mix_ups=CopyNumberReportQC_Task.CopyNumQCMixUps
        # Picard Multiple Metrics Task - NORMAL BAM
        File? normal_bam_picard=normalMM_Task.picard_files
        File? tumor_bam_picard=tumorMM_Task.picard_files
        File? normal_bam_hsmetrics=normalMM_Task.hsMetrics
        File? tumor_bam_hsmetrics=tumorMM_Task.hsMetrics
        File? tumor_cleaned_unmapped_bam=tumorMM_Task.bam_unmapped_cleaned
        File? normal_cleaned_unmapped_bam=normalMM_Task.bam_unmapped_cleaned
        # Cross-Sample Contamination Task
        File contamination_table=select_first([ContEST_Task.contamTable, "null"])
        File normalTable=select_first([ContEST_Task.normTable, "null"])
        File tumorTable=select_first([ContEST_Task.tumTable, "null"])
        Float fracContam=select_first([ContEST_Task.fracContam, "0"])
        # Cross Check Lane Fingerprints Task
        File? cross_check_fingprt_metrics=CrossCheckLaneFingerprints_Task.crossCheckMetrics
    }
}

task CopyNumberReportQC_Task {
    input{
        # TASK INPUT PARAMS
        File tumorBam
        File tumorBamIdx
        File? normalBam
        File? normalBamIdx
        File regionFile
        File readGroupBlackList    
        File captureNormalsDBRCLZip

        # FILE SIZE
        Int tumorBam_size
        Int normalBam_size

        # RUNTIME INPUT PARAMS
        String preemptible = "1"
        String diskGB_boot = "15"
        String diskGB_buffer = "20"
        String memoryGB = "7"
        String cpu = "1"
    }

    # COMPUTE MEMORY SIZE

    Int command_memoryGB = floor(memoryGB) - 1
    
    # COMPUTE DISK SIZE
    Int diskGB = ceil(2.5*ceil(tumorBam_size + normalBam_size + size(regionFile, 'G') + size(readGroupBlackList, 'G') 
                    + size(captureNormalsDBRCLZip, 'G') + diskGB_buffer))

    parameter_meta {
        tumorBam : "sample tumor  BAM file"
        tumorBamIdx : "sample tumor BAI file (BAM indexed)"
        regionFile : ""
        readGroupBlackList : ""
        captureNormalsDBRCLZip : ""
    }

    command {

        # e - exit when a command fails
        # u - exit when script tries to use undeclared variables
        # x - trace what gets executed
        # o - exit status of the last command that threw a non-zero exit code is returned
        set -eux pipefail

        # Copy Number QC for Capture
        # Make lane lists for tumor and normal
        # MakeLaneList_12        
        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/MakeLaneList.jar ${tumorBam}  case.lanelist.txt ;
        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/MakeLaneList.jar ${normalBam} control.lanelist.txt ;

        # make region coverages per lane for tumor
        # RegionCovPerLane_18
        for CHROM in `seq 1 24` ; do
            echo "Working on $CHROM for ${tumorBam}"
            SCATTER_DIR="scatter.$CHROM" ;
            mkdir -v $SCATTER_DIR
            OUTPATH=`echo -ne "$SCATTER_DIR/chr$CHROM.case.rcl"` ;
            echo "OUTPATH is $OUTPATH"
            echo "chrom is $CHROM"
            java -jar /usr/local/bin/RegionCovPerLane.jar ${tumorBam} ${regionFile} $OUTPATH $CHROM
        done ;
        wait ;
        /usr/local/bin/rcl.gather.sh  case.rcl ;
        rm -rf scatter.* ;

        # make region coverages per lane for control
        # RegionCovPerLane_18
        for CHROM in `seq 1 24` ; do
            echo "Working on $CHROM for ${normalBam}"
            SCATTER_DIR="scatter.$CHROM" ;
            mkdir -v $SCATTER_DIR
            OUTPATH=`echo -ne "$SCATTER_DIR/chr$CHROM.control.rcl"` ;
            echo "OUTPATH is $OUTPATH"
            echo "chrom is $CHROM"
            java -jar /usr/local/bin/RegionCovPerLane.jar ${normalBam} ${regionFile} $OUTPATH $CHROM
        done ;
        wait ;
        /usr/local/bin/rcl.gather.sh  control.rcl ;
        rm -rf scatter.* ;

        # run the Matlab based report
        # CopyNumberQCReport_27
        cp -vfr /CopyNumberQCReport_27/unzip/* .

        # Make a file of files paths
        # This command aims to list the zip files contents, filtering out all but the file paths and then write the paths to a file
python <<CODE
from zipfile import ZipFile

with open('capture_normals_db_wdl', 'w') as writer:
    with ZipFile("${captureNormalsDBRCLZip}", 'r') as f:
        names = f.namelist()
        for name in names:
            writer.write(name + '\n')
CODE

        unzip ${captureNormalsDBRCLZip}
        ./run_fh_CopyNumberQCReport.sh /opt/MATLAB/MATLAB_Compiler_Runtime/v710/ PairCopyNumQCReport \
        case.rcl control.rcl case.lanelist.txt control.lanelist.txt \
        ${readGroupBlackList} ${regionFile} capture_normals_db_wdl NA NA

    }

    runtime {
        docker         : "gcr.io/broad-getzlab-workflows/cga_production_pipeline:v0.2"
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible 
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
        #lane lists
        File tumorBamLaneList="case.lanelist.txt"
        File normalBamLaneList="control.lanelist.txt"
        File tumorRCL="case.rcl"
        File normalRCL="control.rcl"
        #CopyNumQC
        File CopyNumQCReport="report.html"
        File CopyNumQCReportPNG="PairCopyNumQCReport_CopyNumberQC.png"
        Int CopyNumQCMixUps=read_int("num.mixups.txt")
    }
}


task PicardMultipleMetrics_Task {
    input{
    # TASK INPUT PARAMS
    File? bam
    File? bamIndex
    String? sampleName
    File targetIntervals
    File baitIntervals
    File refFasta
    File DB_SNP_VCF
    File DB_SNP_VCF_IDX

    String validationStringencyLevel ="LENIENT"
    String run_clean_sam = false
    # FILE SIZE
    Int bam_size
    Int refFasta_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot ="15"
    String diskGB_buffer = "20"
    String memoryGB ="16"
    String cpu ="4"
    String gatk_docker
    String targetedRun = true
    }
    # DEFAULT VALUES

    # COMPUTE MEMORY SIZE
    Int command_memoryGB = floor(memoryGB) - 1

    # COMPUTE DISK SIZE
    Int diskGB = ceil(bam_size + refFasta_size + db_snp_vcf_size + diskGB_buffer)
  
    parameter_meta {
        bam : "sample (normal or tumor) BAM file"
        bamIndex : "sample (normal or tumor) BAI file (BAM indexed)"
        sampleName : "sample (normal or tumor) name, prefix for output"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        DB_SNP_VCF_IDX : "dbSNP indexed file"
    }

    command {

        if [ "${run_clean_sam}" = true ] ;
            then
            gatk --java-options "-Xmx~{command_memoryGB}g" ValidateSamFile \
            --INPUT ${bam} \
            --OUTPUT "${sampleName}.bam_validation" \
            --MODE VERBOSE \
            --IGNORE_WARNINGS true \
            --REFERENCE_SEQUENCE ${refFasta} \
            --VALIDATION_STRINGENCY ${validationStringencyLevel}

            # Run bam through CleanSam to set MAPQ of unmapped reads to zero
            gatk --java-options "-Xmx~{command_memoryGB}g" CleanSam \
            --INPUT ${bam} \
            --OUTPUT ${sampleName}.unmapped_reads_cleaned.bam
        fi

        gatk --java-options "-Xmx~{command_memoryGB}g" CollectMultipleMetrics \
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
        gatk --java-options "-Xmx~{command_memoryGB}g" ConvertSequencingArtifactToOxoG \
        --INPUT_BASE "${sampleName}.multiple_metrics" \
        --OUTPUT_BASE "${sampleName}.multiple_metrics.converted" \
        --REFERENCE_SEQUENCE ${refFasta} \
        --VALIDATION_STRINGENCY ${validationStringencyLevel}

        mkdir ${sampleName}.Picard_Multiple_Metrics
        mv ${sampleName}.multiple_metrics.* ${sampleName}.Picard_Multiple_Metrics
        #zip up reports for QC Nozzle report
        tar -czvf ${sampleName}.picard_multiple_metrics.tar.gz ${sampleName}.Picard_Multiple_Metrics/

        # Collect WES HS metrics
    if [ "${targetedRun}" = true ];
    then
        gatk --java-options "-Xmx~{command_memoryGB}g" CollectHsMetrics \
        --INPUT ${bam} \
        --BAIT_INTERVALS ${targetIntervals} \
        --TARGET_INTERVALS ${baitIntervals} \
        --OUTPUT "${sampleName}.HSMetrics.txt" \
        --VALIDATION_STRINGENCY ${validationStringencyLevel}
    fi 

    }

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible 
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
        File picard_files="${sampleName}.picard_multiple_metrics.tar.gz"
        File? bam_validation="${sampleName}.bam_validation"
        File? bam_unmapped_cleaned="${sampleName}.unmapped_reads_cleaned.bam"
        #File metricsReportsZip="${sampleName}.picard_multiple_metrics.zip"
        #File alignment_summary_metrics="${sampleName}.multiple_metrics.alignment_summary_metrics"
        #File bait_bias_detail_metrics="${sampleName}.multiple_metrics.bait_bias_detail_metrics"
        #File bait_bias_summary_metrics="${sampleName}.multiple_metrics.bait_bias_summary_metrics"
        #File base_distribution_by_cycle="${sampleName}.multiple_metrics.base_distribution_by_cycle.pdf"
        #File base_distribution_by_cycle_metrics="${sampleName}.multiple_metrics.base_distribution_by_cycle_metrics"
        #File gc_bias_detail_metrics="${sampleName}.multiple_metrics.gc_bias.detail_metrics"
        #File gc_bias="${sampleName}.multiple_metrics.gc_bias.pdf"
        #File gc_bias_summary_metrics="${sampleName}.multiple_metrics.gc_bias.summary_metrics"
        #File insert_size_histogram="${sampleName}.multiple_metrics.insert_size_histogram.pdf"
        #File insert_size_metrics="${sampleName}.multiple_metrics.insert_size_metrics"        
        #File pre_adapter_detail_metrics="${sampleName}.multiple_metrics.pre_adapter_detail_metrics"
        #File pre_adapter_summary_metrics="${sampleName}.multiple_metrics.pre_adapter_summary_metrics"
        #File quality_by_cycle="${sampleName}.multiple_metrics.quality_by_cycle.pdf"
        #File quality_by_cycle_metrics="${sampleName}.multiple_metrics.quality_by_cycle_metrics"
        #File quality_distribution="${sampleName}.multiple_metrics.quality_distribution.pdf"
        #File quality_distribution_metrics="${sampleName}.multiple_metrics.quality_distribution_metrics"
        #File quality_yield_metrics="${sampleName}.multiple_metrics.quality_yield_metrics"
        #File converted_oxog_metrics="${sampleName}.multiple_metrics.converted.oxog_metrics"
        File? hsMetrics="${sampleName}.HSMetrics.txt"
    }
}


task ContEST_Task {
    input{    
    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File? normalBam
    File? normalBamIdx
    File targetIntervals
    File gnomad_idx
    File gnomad
    String pairName
    String gatk_docker

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot = "15"
    String diskGB_buffer ="20"
    String memoryGB ="10"
    String cpu ="1"
}
    # DEFAULT VALUES
    # COMPUTE MEMORY SIZE
    Int command_memoryGB = floor(memoryGB) - 1

    # COMPUTE DISK SIZE
    Int diskGB = ceil(tumorBam_size + normalBam_size 
                + size(targetIntervals, "G") + size(gnomad, "G") + diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "sample name, prefix for output"
        targetIntervals : ""
    }

    command <<<

        set -euxo pipefail
        savTable="~{pairName}.contamination.table"

        gatk --java-options "-Xmx~{command_memoryGB}g" GetPileupSummaries \
        -I ~{normalBam} \
        -V ~{gnomad} \
        -L ~{targetIntervals} \
        -O normal.pileups.table

        gatk --java-options "-Xmx~{command_memoryGB}g" GetPileupSummaries \
        -I ~{tumorBam} \
        -V ~{gnomad} \
        -L ~{targetIntervals} \
        -O tumor.pileups.table

        gatk --java-options "-Xmx~{command_memoryGB}g" CalculateContamination \
        -I tumor.pileups.table \
        -matched normal.pileups.table \
        -O $savTable

        awk 'NR == 2 {print $2}' $savTable > contamination.value


    >>>

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible 
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
        File contamTable="${pairName}.contamination.table"
        File normTable="normal.pileups.table"
        File tumTable="tumor.pileups.table"
        Float fracContam=read_float("contamination.value")

    }
}


task CrossCheckLaneFingerprints_Task {
    ## NOTE: UPDATE THE GATK_JAR HERE!
    ## NOTE: Check if there's an update haplotypeDBforcrosscheck
    ## NOTE: Need to update
    input {
    # TASK INPUT PARAMS
    File tumorBam
    File? normalBam
    File tumorBamIdx
    File? normalBamIdx
    String pairName
    File HaplotypeDBForCrossCheck
    String validationStringencyLevel = "LENIENT"
    String gatk_docker

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot ="15"
    String diskGB_buffer ="20"
    String memoryGB ="3"
    String cpu ="1"

	}
 
    # COMPUTE MEMORY SIZE

    Int command_memoryGB = floor(memoryGB) - 1

    # COMPUTE DISK SIZE
    Int diskGB = ceil(tumorBam_size + normalBam_size + size(HaplotypeDBForCrossCheck, "G")
                    + diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        HaplotypeDBForCrossCheck : ""
        validationStringencyLevel : ""
    }

    command {

        set -euxo pipefail

        #drop from haplotypeDB seq entries which aren't in BAM if there are any found
        #PREPPED_HAPLOTYPE_DB=PreppedHaplotypeDB.txt
        #/usr/local/bin/filter_not_in_bam_dict.pl ${normalBam} ${HaplotypeDBForCrossCheck} $PREPPED_HAPLOTYPE_DB
        #CrosscheckLaneFingerprints[version=9]/PREPPED HP DB
        mkdir -v tmp
        gatk --java-options "-Xmx${command_memoryGB}g" CrosscheckFingerprints \
        -I ${tumorBam} \
        -I ${normalBam} \
        -H ${HaplotypeDBForCrossCheck} \
        --TMP_DIR `pwd`/tmp \
        --QUIET false \
        --EXIT_CODE_WHEN_MISMATCH 0 \
        --OUTPUT crosscheck.metrics \
        --VALIDATION_STRINGENCY ${validationStringencyLevel} 

        #Produce crosscheck.stats.txt file for making the html report
        grep -v "#" crosscheck.metrics | sed 1d > crosscheck.metrics.stripped
        
        # python /usr/local/bin/crosscheck_report.py crosscheck.metrics.stripped

    }

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible 
        cpu            : cpu
        disks          : "local-disk ${diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
        File crossCheckMetrics="crosscheck.metrics"
    }
}
