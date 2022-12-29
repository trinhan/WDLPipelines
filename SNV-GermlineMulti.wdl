version 1.0

## Pipeline for running mutation calling from germline samples
## Steps:
## 0. Optional QC using picard
## 1. Strelka2
## 2. Pisces
## 3. CNN variant filer on haplotype caller 
## 4. VEP

#import "https://raw.githubusercontent.com/gatk-workflows/gatk4-germline-snps-indels/master/haplotypecaller-gvcf-gatk4.wdl" as HaplotypeCaller
import "cnn_variant_wdl/cram2filtered.wdl" as CNNFilter
import "pisces_Germline_parallel.wdl" as pisces
import "VEP104.wdl" as VEP
import "run_QC_checks.wdl" as runQC
#import "oncokb.wdl" as oncokb
import "vardict.wdl" as vardict

workflow runGermlineVariants{
    input {
    File normalBam
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File normalBamIdx
    # a string for the name of the pair under analysis used for naming output files
    ##String pairName
    # a string for the name of the tumor sample under analysis used for naming output files
    String ctrlName
    # list of read groups to exclude from the analysis in MuTect1 and MuTect_FC tasks
    File refFasta
    # the FASTA file index for the reference genome (see http://www.htslib.org/doc/faidx.html)
    File refFastaIdx
    # the FASTA file dictionary for the reference genome (see https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)
    File refFastaDict
    # an interval list file that contains the locations of the targets
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

    ## List which callers to use
    Boolean callVardict 
    Boolean callHaplotype     
    Boolean callStrelka 
    Boolean callPisces 

#    File DB_SNP_VCF
#    File DB_SNP_VCF_IDX

#    Boolean runQC
#    Boolean targetedRun

#    String oncotree
#    File AAlist
#    String? grepRm
#    String OncoKBtoken

    ## jointdiscovery inputs
    Array[File] HC_resources
    Array[File] HC_resources_index
    File? gnomad
    File? gnomadidx
    String info_key = "CNN_1D"
    String tensor_type = "reference"
    String cnn_extra_args = "-stand-call-conf 0 -A Coverage -A ChromosomeCounts -A BaseQuality -A FragmentLength -A MappingQuality -A ReadPosition "
    Int? HC_shard_counts

    Int minCallerSupport = 1
    }

    String assembly = if refGenome=="hg19" then "GRCh37" else "GRCh38"

    String targName=basename(sub(targetIntervals,"\\.interval_list", ""))

    Int normalBam_size  = ceil(size(normalBam,  "G") + size(normalBamIdx,   "G")) 
    Int refFasta_size   = ceil(size(refFasta,   "G") + size(refFastaDict,   "G") + size(refFastaIdx, "G")) 
#    Int db_snp_vcf_size = ceil(size(DB_SNP_VCF, "G")+size(DB_SNP_VCF_IDX, "G"))


#    if (runQC){
#        call runQC.PicardMultipleMetrics_Task as normalMM_Task {
#            input:
#                bam=normalBam,
#                bamIndex=normalBamIdx,
#                sampleName=ctrlName,
#                refFasta=refFasta,
#                DB_SNP_VCF=DB_SNP_VCF,
#                DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
#                targetIntervals=targetIntervals,
#                baitIntervals=targetIntervals,
#                gatk_docker=gatk_docker,
#                refFasta_size=refFasta_size,
#                db_snp_vcf_size=db_snp_vcf_size,
#                bam_size=normalBam_size,
#                targetedRun=targetedRun
#
#        }
#    }

        # PREPARE FOR SCATTER
###if ( callVardict || callPisces ){
    call CallSomaticMutations_Prepare_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            gatk_docker=gatk_docker # takes padded interval file (10bp on each side)
    }
###}

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

    if ( callStrelka ){
    call Strelka2Germline_Task {
        input: 
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            normalBam_size=normalBam_size,
            refFasta_size=refFasta_size,
            callRegionsBED=IntervalToBed.output_bed,
            callRegionsBEDTBI=IntervalToBed.output_bed_tbi,
            name=ctrlName
    }
    }

    if ( callVardict ) {
        call vardict.VardictWF as VardictWF {
            input: 
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            tumorBam=normalBam,
            tumorBamIdx=normalBamIdx,
            ctrlName=ctrlName,
            gatk_docker=gatk_docker,
            scatterIndices_in = CallSomaticMutations_Prepare_Task.scatterIndices,
            bed_list_in = CallSomaticMutations_Prepare_Task.bed_list
        }
    }

    if ( callPisces ) {
        call pisces.pisces_workflow as piscesWF {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            normalBam=normalBam,
            normalBai=normalBamIdx,
            pairName=ctrlName,
            pisces_reference=pisces_reference,
            gatk_docker=gatk_docker,
            scatterIndices_in = CallSomaticMutations_Prepare_Task.scatterIndices,
            bed_list_in = CallSomaticMutations_Prepare_Task.bed_list
    }
    }


#    Array[File] VardictFiles = select_first([runvardict.vcfFile, "NULL"])
#    Array[File] PiscesFiles = select_first([runpisces.normal_variants, "NULL"])

    if ( callHaplotype ){
    call CNNFilter.Cram2FilteredVcf as CNNScoreVariantsWorkflow {
        input: 
            input_file=normalBam,                  # Aligned CRAM file or Aligned BAM files
            input_file_index=normalBamIdx,           # Index for an aligned BAM file if that is the input, unneeded if input is a CRAM
            reference_fasta =refFasta,
            reference_dict=refFastaDict,
            reference_fasta_index=refFastaIdx,
            resources = HC_resources,
            resources_index =HC_resources_index,
            output_prefix = ctrlName,
            info_key=info_key,                 # The score key for the INFO field of the vcf (e.g. CNN_1D, CNN_2D)
            tensor_type=tensor_type,              # What kind of tensors the Neural Net expects (e.g. reference, read_tensor)
            scatter_count =select_first([HC_shard_counts, 4]),              # Number of shards for parallelization of HaplotypeCaller and CNNScoreVariants
            snp_tranches=" --snp-tranche 99.9 ",             # Filtering threshold(s) for SNPs in terms of sensitivity to overlapping known variants in resources
            indel_tranches=" --indel-tranche 99.5 " ,          # Filtering threshold(s) for INDELs in terms of sensitivity to overlapping known variants in resources
            gatk_docker=gatk_docker,
            calling_intervals=targetIntervals,
            extra_args=cnn_extra_args
    }
    }

    call Merge_Variants_Germline {
        input:
            ctrlName=ctrlName,
            Haplotype=CNNScoreVariantsWorkflow.cnn_filtered_vcf,
            STRELKA2=Strelka2Germline_Task.strelka2GermlineVCF,
            PISCES_NORMAL=piscesWF.normal_variants,
            Vardict=VardictWF.vardict,
            minCallers=minCallerSupport,
            callVardict=callVardict,
            callStrelka=callStrelka,
            callPisces=callPisces,
            callHaplotype=callHaplotype
    }

    call VEP.variant_effect_predictor as vep {
        input:
          inputFile =Merge_Variants_Germline.MergedGermlineVcf,
          sample_name = ctrlName,
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
          biotype=true,
          pick=true 
    }

    output {
        # Strelka2Germline
       #File? Picard_QC_Output=normalMM_Task.picard_files
       #File? Picard_HsMetrics=normalMM_Task.hsMetrics
       #File? bam_cleaned = normalMM_Task.bam_unmapped_cleaned
       File? strelka2GermlineVCF=Strelka2Germline_Task.strelka2GermlineVCF
       # pisces outputs
       File? pisces_normal_variants=piscesWF.normal_variants
       File? vardict=VardictWF.vardict
       File? HaplotypeVcf=CNNScoreVariantsWorkflow.cnn_filtered_vcf
       File? HaplotypeVcfTbi=CNNScoreVariantsWorkflow.cnn_filtered_vcf_index
      # merged germline output
       File Merged_germline=Merge_Variants_Germline.MergedGermlineVcf
       File Merged_germlineIdx=Merge_Variants_Germline.MergedGermlineVcfIdx
       ## CNN
       File vep_annot = vep.annotatedFile
       File? vep_summary_html=vep.summary_html
     }
}

task Strelka2Germline_Task{
    input {
    File normalBam
    File normalBamIdx
    File refFasta
    File refFastaIdx
    File? callRegionsBED
    File? callRegionsBEDTBI
    String name
    String tmpDIR = "strelkaTMP_" + name
    
     # FILE SIZE
    Int normalBam_size
    Int refFasta_size
    String defthreads ="4"

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot ="15"
    String diskGB_buffer ="20"
    String machine_memoryGB ="24"
    String cpu ="1"
}
    # DEFAULT VALUES

    Int command_memoryGB = ceil(machine_memoryGB) - 1
   
    # COMPUTE DISK SIZE
    Int diskGB = ceil(normalBam_size + refFasta_size + diskGB_buffer)


    command {
        mkdir ${tmpDIR} && /usr/local/bin/configureStrelkaGermlineWorkflow.py \
            --bam ${normalBam} \
            --referenceFasta ${refFasta} \
            --exome \
            --runDir ${tmpDIR} --runDir ${tmpDIR} ${"--callRegions " + callRegionsBED} && \
            ${tmpDIR}/runWorkflow.py -m local -j ${defthreads} && \
            mv ${tmpDIR}/results/variants/variants.vcf.gz  ${name}.strelka2.germline.vcf.gz && \
            mv ${tmpDIR}/results/variants/genome.S1.vcf.gz ${name}.strelka2.genome.germline.vcf.gz


    }
    
    runtime{
        docker : "erictdawson/strelka2:2021-Jan-12"
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        cpu            : cpu 
        disks          : "local-disk ${diskGB} HDD"
        memory         : machine_memoryGB + "GB"
    }
    
    output{
        File strelka2GermlineVCF = "${name}.strelka2.germline.vcf.gz"
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

task Merge_Variants_Germline {

 input {
    # TASK INPUT PARAMS
    File? PISCES_NORMAL
    File? Haplotype
    File? STRELKA2
    File? Vardict
    String ctrlName
    
    # RUNTIME INPUT PARAMS
    String preemptible = "2"
    String diskGB_boot = "10"
    String diskGB_buffer ="5"
    String memoryGB ="4"
    String cpu ="1"
    Int minCallers

    Boolean callPisces
    Boolean callVardict
    Boolean callStrelka
    Boolean callHaplotype

    }

    # DEFAULT VALUES
    Int minV = minCallers - 1
       
    Int diskGB = ceil(size(Haplotype, "G"))*4  +  diskGB_buffer

    command <<<

        RenameFiles=""

        # reformat if strelka has been run
        if [ ~{callStrelka} == true ]; then
        STRELKA_pass="~{ctrlName}.S2.PASS.vcf"
        bcftools view -f PASS "~{STRELKA2}" > $STRELKA_pass
        bgzip $STRELKA_pass
        tabix -p vcf $STRELKA_pass.gz
        RenameFiles="${RenameFiles}~{ctrlName}.Strelka\n"
        fi

        # reform for piscse
        if [ ~{callPisces} == true ]; then
        PISCES_pass="~{ctrlName}.Pisces.pass.vcf"
        bcftools view -f PASS "~{PISCES_NORMAL}" > $PISCES_pass
        sed -i 's/##FORMAT=<ID=AD,Number=R,/##FORMAT=<ID=AD,Number=.,/g' $PISCES_pass
        bgzip $PISCES_pass
        tabix -p vcf $PISCES_pass.gz
        RenameFiles="${RenameFiles}~{ctrlName}.Pisces\n"
        fi

        # reformat for vardict
        if [ ~{callVardict} == true ]; then
        Vardict_PASSED="~{ctrlName}.Vardict.passed.vcf"
        bcftools view -f PASS "~{Vardict}" > $Vardict_PASSED
        sed -i 's/##FORMAT=<ID=AD,Number=R,/##FORMAT=<ID=AD,Number=.,/g' $Vardict_PASSED
        bgzip $Vardict_PASSED
        tabix -p vcf $Vardict_PASSED.gz
        RenameFiles="${RenameFiles}~{ctrlName}.Vardict\n"
        fi

        # reformat for Haplotupe
        if [ ~{callHaplotype} == true ]; then
        HP_pass="~{ctrlName}.haplo.pass.vcf"
        bcftools view -f PASS "~{Haplotype}" > $HP_pass
        sed -i 's/##FORMAT=<ID=AD,Number=R,/##FORMAT=<ID=AD,Number=.,/g' $HP_pass
        bgzip $HP_pass
        tabix -p vcf $HP_pass.gz
        RenameFiles="${RenameFiles}~{ctrlName}.Haplotype\n"
        fi

        echo -e $RenameFiles > samples.txt

         #STRELKA_unzip="~{ctrlName}.S2.unzip.vcf"
        #HP_unzip="~{ctrlName}.haplo.vcf"

        MERGED_VCF="~{ctrlName}.SNV.mergedGermline.vcf.gz"
        RENAME_MERGED_VCF_ALL="~{ctrlName}.mergedGermline.all.vcf.gz"
        RENAME_MERGED_VCF_ANN="~{ctrlName}.mergedGermline.ann.vcf"
        RENAME_MERGED_VCF_FILT="~{ctrlName}.mergedGermline.filt.vcf"

        #merge vcfs
        bcftools merge ~{true="$STRELKA_pass.gz" false="" callStrelka} \
                       ~{true="$PISCES_pass.gz" false="" callPisces} \
                       ~{true="$Vardict_PASSED.gz" false="" callVardict} \
                       ~{true="$HP_pass.gz" false="" callHaplotype} -O vcf -o $MERGED_VCF --force-samples
        bcftools reheader -s samples.txt $MERGED_VCF > $RENAME_MERGED_VCF_ALL
        tabix -p vcf $RENAME_MERGED_VCF_ALL

        #merge
        bcftools query --format '%CHROM\t%POS\t%POS\n' $RENAME_MERGED_VCF_ALL > test.output 
        bcftools query -f '[\t%SAMPLE=%GT]\n' $RENAME_MERGED_VCF_ALL | awk '{print 4-gsub(/.\/\./, "")}' > output
        paste test.output output > annots.tab 
        bgzip annots.tab
        tabix -s1 -b2 -e2 annots.tab.gz
        echo '##INFO=<ID=NCALLS,Number=1,Type=Integer,Description="Number of callers">' > annots.hdr
        bcftools annotate -a annots.tab.gz -h annots.hdr -c CHROM,FROM,TO,NCALLS $RENAME_MERGED_VCF_ALL > $RENAME_MERGED_VCF_ANN

        bgzip $RENAME_MERGED_VCF_ANN
        tabix -p vcf $RENAME_MERGED_VCF_ANN.gz
        # querybased on the number of callers supporting
        echo 'filter ino'
        bcftools filter -i'NCALLS>~{minV}' $RENAME_MERGED_VCF_ANN.gz -o $RENAME_MERGED_VCF_FILT
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
        File MergedGermlineVcf="~{ctrlName}.mergedGermline.filt.vcf.gz"
        File MergedGermlineVcfIdx="~{ctrlName}.mergedGermline.filt.vcf.gz.tbi"
    }     
 
}

task CombineVariants {
    input {
        Array[File] input_header
        String caller
        File ref_fasta
        File ref_fai
        File ref_dict
        # runtime
        String gatk_docker
        String sample_name
        Int mem_gb= 6
    }
        Int diskGB = 4*ceil(size(ref_fasta, "GB")+size(input_header, "GB"))

    command <<<
       gatk GatherVcfs -I ~{sep=' -I ' input_header} -R ~{ref_fasta} -O ~{sample_name}.~{caller}.vcf
    >>>

    runtime {
        docker: gatk_docker
        memory: "~{mem_gb} GB"
        disk_space: "local-disk ~{diskGB} HDD"
    }

    output {
    File merged_vcf = "~{sample_name}.~{caller}.vcf"
    }
}

task UpdateHeaders {
    input {
        Array[File] input_vcfs
        File ref_dict
        # runtime
        String gatk_docker
        String caller
        Int mem_gb=6
    }
        Int diskGB = 4*ceil(size(ref_dict, "GB"))

    command <<<

        count=0
        for i in ~{sep=' ' input_vcfs}; 
        do 
         newstr=`basename $i`
         gatk UpdateVCFSequenceDictionary \
            -V $i \
            --source-dictionary ~{ref_dict} \
            --output $newstr.$count.reheader.~{caller}.vcf \
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
    Array[File] head_vcf = glob("*.reheader.~{caller}.vcf")
    }
}


