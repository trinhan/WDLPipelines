version 1.0

# this runs mutect1 mutect2 strelka2

import "mutect2_wdl/mutect2.wdl" as Mutect2WF
import "vardict.wdl" as vardict

workflow runVariantCallers{
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
    # the FASTA file index for the reference genome (see http://www.htslib.org/doc/faidx.html)
    File refFastaIdx
    # the FASTA file dictionary for the reference genome (see https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)
    File refFastaDict
    # an interval list file that contains the locations of the targets
    File targetIntervals
    File? gnomad
    File? gnomad_idx
    File? m2_extra_args
    File? m2_pon
    File? m2_pon_idx
    File? variants_for_contamination
    File? variants_for_contamination_idx
    File? DB_SNP_VCF
    # index file of VCF file of DB SNP variants
    File? DB_SNP_VCF_IDX
    # catalogue of somatic mutations in VCF format    
    File? cosmicVCF # can this be updated?
    # contamination fraction: from running ContEST  
    Float fracContam 
    File strelka_config
    # Options here: Pair TumorOnly, Germline
    String runMode = if defined(normalBam) then "Paired" else "TumOnly"
    String gatk_docker
    Int? minCallerSupport

    }

    String targName=basename(sub(targetIntervals,"\\.interval_list", ""))

    Boolean runS2 = if (runMode =="Paired") then true else false

    Int tumorBam_size   = ceil(size(tumorBam,   "G") + size(tumorBamIdx,    "G")) 
    Int normalBam_size  = if defined (normalBam) then ceil(size(normalBam,  "G") + size(normalBamIdx,   "G")) else 0
    Int db_snp_vcf_size = if defined (DB_SNP_VCF) then ceil(size(DB_SNP_VCF, "G") + size(DB_SNP_VCF_IDX, "G")) else 0
    Int refFasta_size   = ceil(size(refFasta,   "G") + size(refFastaDict,   "G") + size(refFastaIdx, "G")) 



        # PREPARE FOR SCATTER
    call CallSomaticMutations_Prepare_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            gatk_docker=gatk_docker # takes padded interval file (10bp on each side)
    }

    call CreateFoFN {
        input:
            array_of_files = CallSomaticMutations_Prepare_Task.interval_files,
            fofn_name = targName,
            docker = gatk_docker
    }

    call IntervalToBed {
        input:
            targetIntervals=targetIntervals,
            gatk_docker=gatk_docker,
            output_name = targName
    }

    #SCATTER AND ANALYZE
    scatter (idx in CallSomaticMutations_Prepare_Task.scatterIndices) {
        # Identification of somatic point mutations in next generation sequencing data of cancer genomes.
        call Mutect1_Task {
            input:
                tumorBam=tumorBam,
                tumorBamIdx=tumorBamIdx,
                normalBam=normalBam,
                normalBamIdx=normalBamIdx,
                pairName=pairName,
                caseName=caseName,
                ctrlName=ctrlName,
                fracContam=fracContam,
                mutectIntervals=CallSomaticMutations_Prepare_Task.interval_files[idx],
                refFasta=refFasta,
                refFastaIdx=refFastaIdx,
                refFastaDict=refFastaDict,
                DB_SNP_VCF=DB_SNP_VCF,
                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
                cosmicVCF=cosmicVCF,
                MuTectNormalPanel=m2_pon,
                MuTectNormalPanelIdx=m2_pon_idx,
                refFasta_size=refFasta_size,
                db_snp_vcf_size=db_snp_vcf_size,
                tumorBam_size=tumorBam_size,
                normalBam_size=normalBam_size
        }

        call vardict.VarDict as runvardict {
            input:
                referenceFasta=refFasta,
                referenceFastaFai=refFastaIdx,
                tumorBam=tumorBam,
                normalBam=normalBam,
                normalBamIndex=normalBamIdx,
                tumorBamIndex=tumorBamIdx,
                outputName=pairName,
                bedFile=CallSomaticMutations_Prepare_Task.bed_list[idx],
                tumorSampleName=caseName,
                normalSampleName=ctrlName
        }
    }

    call Mutect2WF.Mutect2 as M2WF {
        input:
            intervals=targetIntervals,
            ref_fasta=refFasta,
            ref_fai=refFastaIdx,
            ref_dict=refFastaDict,
            tumor_reads=tumorBam,
            tumor_reads_index=tumorBamIdx,
            normal_reads=normalBam,
            normal_reads_index=normalBamIdx,
            pon=m2_pon,
            pon_idx=m2_pon_idx,
            scatter_count=24,
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            variants_for_contamination=variants_for_contamination,
            variants_for_contamination_idx=variants_for_contamination_idx,
            m2_extra_args=m2_extra_args,
            gatk_docker=gatk_docker
        }

    call UpdateHeaders {
        input:
            input_vcfsVar =runvardict.vcfFile,
            ref_dict=refFastaDict,
            gatk_docker = gatk_docker
    }

    call CombineVariants {
        input:
            input_VD=UpdateHeaders.VDhead_vcf,
            ref_fasta = refFasta,
            ref_fai = refFastaIdx,
            ref_dict = refFastaDict,
            gatk_docker = gatk_docker,
            sample_name = caseName
    }    

    if (runS2){
    call Strelka2Somatic_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            callRegionsBED=IntervalToBed.output_bed,
            callRegionsBEDTBI=IntervalToBed.output_bed_tbi,
            tumorBam_size=tumorBam_size,
            normalBam_size=normalBam_size,
            refFasta_size=refFasta_size,
            name=pairName,
            config=strelka_config
    }
    }


    if (runMode=="Paired") {
        call Merge_Variant_Calls as PairedCall {
            input:
            pairName=pairName,
            mutect1_cs=Mutect1_Task.mutect1_cs,
            M2=M2WF.filtered_vcf,
            STRELKA2_SNVS=Strelka2Somatic_Task.strelka2SomaticSNVs,
            STRELKA2_INDELS=Strelka2Somatic_Task.strelka2SomaticIndels,
            Vardict=CombineVariants.merged_vcfVD,
            ctrlName=ctrlName,
            caseName=caseName,
            refFastaDict=refFastaDict,
            runMode=runMode
        } 
    } 

    if (runMode=="TumOnly") {
        call Merge_Variant_Calls as TumCall {
            input:
            pairName=pairName,
            mutect1_cs=Mutect1_Task.mutect1_cs,
            M2=M2WF.filtered_vcf,
            caseName=caseName,
            Vardict=CombineVariants.merged_vcfVD,
            refFastaDict=refFastaDict,
            runMode=runMode
        } 
    }

     output {
        
       File? strelka2SomaticSNVs = Strelka2Somatic_Task.strelka2SomaticSNVs
       File? strelka2SomaticIndels = Strelka2Somatic_Task.strelka2SomaticIndels
        File vardict_out=CombineVariants.merged_vcfVD
      #  M2 workflow2 outputs
       File M2_filtered_vcf=M2WF.filtered_vcf
       File M2_filtered_vcf_idx=M2WF.filtered_vcf_idx
      #  merged output haplotypecaller
       File Combined_raw_variants_gz=select_first([PairedCall.MergedVcfGz, TumCall.MergedVcfGz])
       File Combined_raw_variants_tbi=select_first([PairedCall.MergedVcfIdx, TumCall.MergedVcfIdx])
       File Combined_raw_variants_maf=select_first([PairedCall.MergedMaf, TumCall.MergedMaf])
       File VariantSitesBed=select_first([PairedCall.LocBed, TumCall.LocBed])
        }
}

task Mutect1_Task {
    input{
    # TASK INPUT PARAMS
    File tumorBam
    File? normalBam
    File tumorBamIdx
    File? normalBamIdx
    String pairName
    String caseName
    String? ctrlName
    File mutectIntervals
    File? DB_SNP_VCF
    File? DB_SNP_VCF_IDX
    File? cosmicVCF
    File? MuTectNormalPanel
    File? MuTectNormalPanelIdx
    File refFasta
    File refFastaIdx
    File refFastaDict
    Float fracContam

    String downsample ="9999"

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size
    Int db_snp_vcf_size

    # RUNTIME INPUT PARAMS
    String preemptible = "1"
    String diskGB_boot = "15"
    String diskGB_buffer = "20"
    String machine_memoryGB ="15"
    String cpu ="1"
}

    # COMPUTE MEMORY SIZE

    Int command_memoryGB = ceil(machine_memoryGB) - 1
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + db_snp_vcf_size
                + size(mutectIntervals, "G") + size(cosmicVCF, "G") + size(MuTectNormalPanel, "G") + diskGB_buffer)
    

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        pairName : "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
        mutectIntervals : "a list of genomic intervals over which MuTect1 will operate"
        DB_SNP_VCF : "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis by some PROGRAMs"
        cosmicVCF : "catalogue of somatic mutations in VCF format"
        MuTectNormalPanel : "1000 genomes panel of normals in VCF format"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
        fracContam : "fraction of cross-sample contamination, output from ContEst task"
        downsample : "downsample reads to a given capping threshold coverage"
    }

    command {

        set -euxo pipefail

        #variable for normal panel

        java "-Xmx${command_memoryGB}g" -jar /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
        -L ${mutectIntervals} \
        ${"--normal_sample_name " + ctrlName} \
        ${"-I:normal " + normalBam} \
        --tumor_sample_name ${caseName} \
        -I:tumor ${tumorBam} \
        --reference_sequence ${refFasta} \
        --fraction_contamination ${fracContam} \
        ${"--dbsnp " + DB_SNP_VCF} \
        ${"--cosmic " + cosmicVCF} \
        ${"--normal_panel " + MuTectNormalPanel} \
        --out ${pairName}.MuTect1.call_stats.txt \
        --coverage_file ${pairName}.MuTect1.coverage.wig.txt \
        --power_file ${pairName}.MuTect1.power.wig.txt \
        --downsample_to_coverage ${downsample} 

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
        File mutect1_cs="${pairName}.MuTect1.call_stats.txt"
        File mutect1_pw="${pairName}.MuTect1.power.wig.txt"
        File mutect1_cw="${pairName}.MuTect1.coverage.wig.txt"
    }
}

task Strelka2Somatic_Task {
    input {
    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File? normalBam
    File? normalBamIdx
    File refFasta
    File refFastaIdx
    File config
    File? callRegionsBED
    File? callRegionsBEDTBI
    String name
   
    String tmpDIR = "strelkaTMP_" + name

    # FILE SIZE
    Int tumorBam_size
    Int normalBam_size
    Int refFasta_size
    String defthreads ="4"

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot = "15"
    String diskGB_buffer ="20"
    String machine_memoryGB ="25"
    String cpu ="1"
}
    # DEFAULT VALUES
    Int diskGB = ceil(tumorBam_size + normalBam_size + refFasta_size + diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM file)"
        normalBam : "sample normal BAM file"
        normalBamIdx : "sample normal BAI file (indexed BAM file)"
        name : "a string for the name of the pair under analysis used for naming output files"
        refFasta : "FASTA file for reference genome"
        refFastaIdx : "FASTA file index for the reference genome"
        config : "Strelka configuration file"
    }

    command{

        set -euxo pipefail
        
        mkdir ${tmpDIR} && /usr/local/bin/configureStrelkaSomaticWorkflow.py \
            --normalBam ${normalBam} \
            --tumorBam ${tumorBam} \
            --referenceFasta ${refFasta} \
            --exome \
            --runDir ${tmpDIR} ${"--callRegions " + callRegionsBED} && \
            ${tmpDIR}/runWorkflow.py -m local -j ${defthreads} && \
            gunzip -c ${tmpDIR}/results/variants/somatic.snvs.vcf.gz > ${name}.strelka2.somatic.snvs.vcf && \
            gunzip -c ${tmpDIR}/results/variants/somatic.indels.vcf.gz > ${name}.strelka2.somatic.indels.vcf
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
        File strelka2SomaticSNVs = "${name}.strelka2.somatic.snvs.vcf"
        File strelka2SomaticIndels = "${name}.strelka2.somatic.indels.vcf"
    }
}

task CallSomaticMutations_Prepare_Task {
    input {
    # TASK INPUT PARAMS
    File targetIntervals
    File refFasta
    File refFastaIdx
    File refFastaDict
    String nWay = "10"
    # RUNTIME INPUT PARAMS
    String preemptible = "1"
    String diskGB_boot = "15"

    String gatk_docker
}

    parameter_meta {
        nWay : "Number of ways to scatter (MuTect1 and MuTect2)"
        targetIntervals : "a list of genomic intervals over which MuTect1 will operate"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
    }

    command {
        set -euxo pipefail
        seq 0 $((~{nWay}-1)) > indices.dat
        # create a list of intervalfiles
        mkdir intervalfolder
        gatk SplitIntervals -R ~{refFasta} -L ~{targetIntervals} --scatter-count ~{nWay} --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION -O intervalfolder
        cp intervalfolder/*.interval_list .

        ## make the list of bed files
        mkdir bedfolder
        for file in *.interval_list;
        do 
            gatk IntervalListToBed -I $file -O bedfolder/$file.bed
            ##small hack to subtract 1 from the bed file
        done 

        cp bedfolder/*.bed .
    }

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        memory         : "1 GB"
    }

    output {
        Array[File] interval_files=glob("*.interval_list")
        Array[Int] scatterIndices=read_lines("indices.dat")
        Array[File] bed_list=glob("*.bed")
    }
}

task CreateFoFN {
  input {
    # Command parameters
    Array[String] array_of_files
    String fofn_name
    # Runtime parameters
    String docker
  }


  command {
    mv ${write_lines(array_of_files)} "${fofn_name}.interval.list"
  }
  output {
    File fofn_list = "${fofn_name}.interval.list"
  }
  runtime {
    docker: docker
    preemptible: 3
  }
}

task IntervalToBed {
    input {
    File targetIntervals
    String? bed_intervallist_extra_args
    String output_name

    # runtime
    String gatk_docker
    String disk_space = "15"
    String preempt = "2"
    String machine_mem = "4"
    }

    Int command_mem = floor(machine_mem) - 1
    String output_name_intervals = output_name + ".bed"

    command {
        set -e
         
        gatk --java-options "-Xmx${command_mem}g" IntervalListToBed \
            -I ${targetIntervals} \
            -O ${output_name_intervals} ${bed_intervallist_extra_args}

        bgzip ${output_name_intervals}

        tabix -s 1 -b 2 -e 3 "${output_name_intervals}.gz"

        # also set a command here to split the intervals? or grep bed according to chromosome?
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preempt
    }

    output {
        File output_bed = "${output_name_intervals}.gz"
        File output_bed_tbi = "${output_name_intervals}.gz.tbi"
    }
}

task Merge_Variant_Calls {
    input {
        Array[File] mutect1_cs
        File M2
        File? STRELKA2_INDELS
        File? STRELKA2_SNVS
        File Vardict
        String pairName
        String caseName
        String? ctrlName
        File refFastaDict
        String? preemptible = "2"
        String? diskGB_boot = "10"
        String? diskGB_buffer ="5"
        String? memoryGB ="4"
        String? cpu ="1"
        String runMode 
        Int minCallers =2
    }

    # DEFAULT VALUES
    Int minV = minCallers - 1
       
    Int diskGB = ceil(size(M2, "G"))*4  +  diskGB_buffer
    String runS2 =if defined(STRELKA2_INDELS) then "1" else "0"
    
    parameter_meta {
        mutect1_cs : "list of mutect variants"
        M2 : "list of mutect2 variants"
        STRELKA2_INDELS : "list of strelka2 variants"
        STRELKA2_SNVS: "list of strelka2 snvs"
        pairName: "a string for the name of the pair under analysis used for naming output files"
        caseName : "tumor sample name, prefix for output"
        ctrlName : "normal sample name, prefix for output"
    }

    command <<<
        set -x

python3 <<CODE

mutect1_cs_file_array = '~{sep="," mutect1_cs}'.split(",")
print(mutect1_cs_file_array)
mutect1_cs_list = open('mutect1_cs_list.txt', 'w')

for i in range(len(mutect1_cs_file_array)):
    mutect1_cs_list.write(mutect1_cs_file_array[i] + '\n')
mutect1_cs_list.close()

CODE
        
        # MuTect1 files
        MUTECT1_CS="~{pairName}.MuTect1.call_stats.txt"
        MUTECT1_CS_PASSED="~{pairName}.MuTect1.call_stats.passed.txt"
        MUTECT1_CS_REJECTED="~{pairName}.MuTect1.call_stats.rejected.txt"
        MUTECT1_CS_VCF="~{pairName}.MuTect1.call_stats.vcf"
        # MuTect2 files
        MUTECT2_CS_PASSED="~{pairName}.MuTect2.call_stats.passed.vcf"
        M2temp="~{pairName}.MuTect2.call_stats.vcf"
        Vardict_PASSED="~{pairName}.Vardict.passed.vcf"
        python3 /usr/local/bin/merge_callstats.py "mutect1_cs_list.txt" $MUTECT1_CS
        # Filter MuTect1 mutation calls that passed filter
        python3 /usr/local/bin/filter_passed_mutations.py $MUTECT1_CS $MUTECT1_CS_PASSED $MUTECT1_CS_REJECTED "KEEP"
        # Convert MuTect1 call stats to VCF # can revert to the previous verison of this?
        python3 /usr/local/bin/M1_txt2vcf.py $MUTECT1_CS_PASSED $MUTECT1_CS_VCF "~{caseName}" "~{ctrlName}" -VC "None" --runMode ~{runMode}
        bgzip $MUTECT1_CS_VCF
        tabix -p vcf $MUTECT1_CS_VCF.gz
          
        # Name all the merged files
        MERGED_VCF="~{pairName}.M1_M2_S2_vardict.passed.merged.vcf.gz"
        RENAME_MERGED_VCF="~{pairName}.M1_M2_S2_vardict.passed.merged2.vcf.gz"
        RENAME_MERGED_VCF_decomp="~{pairName}.M1_M2_S2_vardict.passed.merged2.vcf"
        RENAME_MERGED_VCF_ANN="~{pairName}.M1_M2_S2_vardict.merged.ann.vcf"
        RENAME_MERGED_VCF_FILT="~{pairName}.M1_M2_S2_vardict.merged.filt.vcf"

        echo 'merge vardict files'
        bgzip ~{Vardict}
        tabix -p vcf ~{Vardict}.gz

    if [ ~{runMode} == "Paired" ];
        then
        echo -e "~{ctrlName}.M1\n~{caseName}.M1\n~{ctrlName}.M2\n~{caseName}.M2\n~{ctrlName}.S2\n~{caseName}.S2\n~{caseName}.Vardict\n~{ctrlName}.Vardict\n"> samples.txt
        echo -e "~{caseName}.M1\n~{caseName}.M2\n~{caseName}.S2\n~{caseName}.Vardict\n"> samples_tum.txt
        elif [ ~{runMode} == "TumOnly" ] ;
        then
        echo -e "~{caseName}.M1\n~{caseName}.M2\n~{caseName}.Vardict\n"> samples.txt ##~{caseName}.P\n
        cp samples.txt samples_tum.txt
    fi;

        echo 'filter out failed variants'

        bcftools view -f PASS ~{M2} > $MUTECT2_CS_PASSED
        bgzip $MUTECT2_CS_PASSED
        tabix -p vcf $MUTECT2_CS_PASSED.gz
        bcftools view -f PASS ~{Vardict}.gz > $Vardict_PASSED
        bgzip $Vardict_PASSED
        tabix -p vcf $Vardict_PASSED.gz
        cat samples.txt
        # Strelka2 files
    if [ ~{runS2} -eq "1" ];
    then
        STRELKA2_INDEL_REFORMATTED_VCF="~{pairName}.Strelka2.call_stats.indel.re_formatted.vcf"
        STRELKA2_INDEL_PASSED="~{pairName}.Strelka2.call_stats.indel.passed.vcf"
        STRELKA2_SNV_PASSED="~{pairName}.Strelka2.call_stats.snv.passed.vcf"
        STRELKA2_SNV_REFORMATTED_VCF="~{pairName}.Strelka2.call_stats.snv.re_formatted.vcf"
        STRELKA2_MERGE="~{pairName}.Strelka2.call_stats.merged.vcf.gz"

        bcftools view -f PASS ~{STRELKA2_INDELS} > $STRELKA2_INDEL_PASSED
        bcftools view -f PASS ~{STRELKA2_SNVS} > $STRELKA2_SNV_PASSED 

        python3 /usr/local/bin/strelka_allelic_count_snv.py $STRELKA2_SNV_PASSED $STRELKA2_SNV_REFORMATTED_VCF "~{caseName}" "~{ctrlName}" -VC "None"

        python3 /usr/local/bin/strelka_allelic_count_indel.py $STRELKA2_INDEL_PASSED $STRELKA2_INDEL_REFORMATTED_VCF "~{caseName}" "~{ctrlName}" -VC "None"
        
        # compress the outputs 
        bgzip $STRELKA2_INDEL_REFORMATTED_VCF
        bgzip $STRELKA2_SNV_REFORMATTED_VCF
        tabix -p vcf $STRELKA2_INDEL_REFORMATTED_VCF.gz
        tabix -p vcf $STRELKA2_SNV_REFORMATTED_VCF.gz

        bcftools concat -a $STRELKA2_SNV_REFORMATTED_VCF.gz $STRELKA2_INDEL_REFORMATTED_VCF.gz -O vcf -o $STRELKA2_MERGE
        tabix -p vcf $STRELKA2_MERGE

        bcftools merge $MUTECT1_CS_VCF.gz $MUTECT2_CS_PASSED.gz $STRELKA2_MERGE $Vardict_PASSED.gz -O vcf -o $MERGED_VCF --force-samples 
    else 
        bcftools merge $MUTECT1_CS_VCF.gz $MUTECT2_CS_PASSED.gz $Vardict_PASSED.gz -O vcf -o $MERGED_VCF --force-samples 
    fi;

        bcftools reheader -s samples.txt $MERGED_VCF > $RENAME_MERGED_VCF
        tabix -p vcf $RENAME_MERGED_VCF

        echo 'filter based on the number of callers'
        bcftools query --format '%CHROM\t%POS\t%POS\n' $RENAME_MERGED_VCF > test.output 
        bcftools query -f '[\t%AF]\n' $RENAME_MERGED_VCF -S samples_tum.txt | awk '{print gsub(/0\./, "")}' > output
        paste test.output output > annots.tab 
        bgzip annots.tab
        tabix -s1 -b2 -e2 annots.tab.gz
        echo '##INFO=<ID=NCALLS,Number=1,Type=Integer,Description="Number of callers">' > annots.hdr
        bcftools annotate -a annots.tab.gz -h annots.hdr -c CHROM,FROM,TO,NCALLS $RENAME_MERGED_VCF > $RENAME_MERGED_VCF_ANN
        bcftools filter -i'NCALLS>~{minV}' $RENAME_MERGED_VCF_ANN -o $RENAME_MERGED_VCF_FILT


        # extract the variant locations for mutect2
        ##gunzip -c $RENAME_MERGED_VCF > $RENAME_MERGED_VCF_decomp
        python3 /usr/local/bin/vcf2mafbed.py $RENAME_MERGED_VCF_FILT "~{pairName}.M1_M2_S2_vardict.passed.filt.maf" "~{pairName}.intervals.bed" 150 "~{runMode}"
        # run picard to change the input 
        java -jar /tmp/picard.jar BedToIntervalList -I "~{pairName}.intervals.bed" -O "~{pairName}.variantList.interval_list" -SD ~{refFastaDict}

        # compress the VCF file as last step
        bgzip $RENAME_MERGED_VCF_FILT
        tabix -p vcf $RENAME_MERGED_VCF_FILT.gz
        >>>
        
    runtime {
        docker         : "trinhanne/sambcfhts:v1.13.3"
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ~{diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
        File MergedVcfGz="~{pairName}.M1_M2_S2_vardict.merged.filt.vcf.gz"
        File MergedVcfIdx="~{pairName}.M1_M2_S2_vardict.merged.filt.vcf.gz.tbi"
        File MergedMaf="~{pairName}.M1_M2_S2_vardict.passed.filt.maf"
        File LocBed="~{pairName}.intervals.bed"
        File Interval_list="~{pairName}.variantList.interval_list"
    }
}

task CombineVariants {
    input {
        Array[File] input_VD
        File ref_fasta
        File ref_fai
        File ref_dict
        # runtime
        String gatk_docker
        String sample_name
        Int mem_gb= 6
    }
        Int diskGB = 4*ceil(size(ref_fasta, "GB")+size(input_VD, "GB"))

    command <<<

        gatk GatherVcfs -I ~{sep=' -I ' input_VD} -R ~{ref_fasta} -O ~{sample_name}.VD.vcf

    >>>

    runtime {
        docker: gatk_docker
        memory: "~{mem_gb} GB"
        disk_space: "local-disk ~{diskGB} HDD"
    }

    output {
    File merged_vcfVD = "~{sample_name}.VD.vcf"
    }
}

task UpdateHeaders {
    input {
        Array[File] input_vcfsVar
        File ref_dict
        # runtime
        String gatk_docker
        Int mem_gb= 6
    }
        Int diskGB = 4*ceil(size(ref_dict, "GB")+size(input_vcfsVar, "GB"))

    command <<<

        # vardict
        count=0
        for i in ~{sep=' ' input_vcfsVar}; 
        do 
         newstr=`basename $i`
         gatk UpdateVCFSequenceDictionary \
            -V $i \
            --source-dictionary ~{ref_dict} \
            --output $newstr.$count.reheaderVD.vcf \
            --replace true
         count+=1
        done

    >>>

    runtime {
        docker: gatk_docker
        memory: "~{mem_gb} GB"
        disk_space: "local-disk ~{diskGB} HDD"
    }

    output {
    Array[File] VDhead_vcf = glob("*.reheaderVD.vcf")
    }
}

