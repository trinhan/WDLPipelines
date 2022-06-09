version 1.0

## Pipeline for running mutation calling from germline samples
## Steps:
## 1. Strelka2
## 2. Pisces
## 3. CNN variant filer on haplotype caller 
## 4. VEP

#import "https://raw.githubusercontent.com/gatk-workflows/gatk4-germline-snps-indels/master/haplotypecaller-gvcf-gatk4.wdl" as HaplotypeCaller
import "cnn_variant_wdl/cram2filtered.wdl" as CNNFilter
import "pisces_task.wdl" as pisces
import "VEP104.wdl" as VEP
import "run_QC_checks.wdl" as runQC
import "oncokb.wdl" as oncokb

workflow runGermlineVariants{
    input {
    File normalBam
    # sample normal BAI file (BAM indexed) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File normalBamIdx
    # a string for the name of the pair under analysis used for naming output files
    String pairName
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

    File DB_SNP_VCF
    File DB_SNP_VCF_IDX

    Boolean runQC
    Boolean targetedRun

    String oncotree
    File AAlist
    String? canonical
    String? grepRm
    String OncoKBtoken

    ## jointdiscovery inputs
    Array[File] HC_resources
    Array[File] HC_resources_index
    File gnomad
    File gnomadidx
    String info_key
    String tensor_type = "reference"
    String cnn_extra_args = "-stand-call-conf 0 -A Coverage -A ChromosomeCounts -A BaseQuality -A FragmentLength -A MappingQuality -A ReadPosition "
    Int? HC_shard_counts
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
            name=pairName
    }

    call pisces.runpisces as runpisces {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaIdx,
            refFastaDict=refFastaDict,
            normalBam=normalBam,
            normalBai=normalBamIdx,
            pairName=pairName,
            pisces_reference=pisces_reference,
            interval=targetIntervals,
            runMode="Germline"   
    }

    call CNNFilter.Cram2FilteredVcf as CNNScoreVariantsWorkflow {
        input: 
            input_file=normalBam,                  # Aligned CRAM file or Aligned BAM files
            input_file_index=normalBamIdx,           # Index for an aligned BAM file if that is the input, unneeded if input is a CRAM
            reference_fasta =refFasta,
            reference_dict=refFastaDict,
            reference_fasta_index=refFastaIdx,
            resources = HC_resources,
            resources_index =HC_resources_index,
            output_prefix = pairName,
            info_key=info_key,                 # The score key for the INFO field of the vcf (e.g. CNN_1D, CNN_2D)
            tensor_type=tensor_type,              # What kind of tensors the Neural Net expects (e.g. reference, read_tensor)
            scatter_count =select_first([HC_shard_counts, 4]),              # Number of shards for parallelization of HaplotypeCaller and CNNScoreVariants
            snp_tranches=" --snp-tranche 99.9 ",             # Filtering threshold(s) for SNPs in terms of sensitivity to overlapping known variants in resources
            indel_tranches=" --indel-tranche 99.5 " ,          # Filtering threshold(s) for INDELs in terms of sensitivity to overlapping known variants in resources
            gatk_docker=gatk_docker,
            calling_intervals=targetIntervals,
            extra_args=cnn_extra_args
    }

    call Merge_Variants_Germline {
        input:
            ctrlName=ctrlName,
            Haplotype=CNNScoreVariantsWorkflow.cnn_filtered_vcf,
            STRELKA2=Strelka2Germline_Task.strelka2GermlineVCF,
            PISCES_NORMAL=runpisces.normal_variants      
    }

    call VEP.variant_effect_predictor as vep {
        input:
          inputFile =Merge_Variants_Germline.MergedGermlineVcf,
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
          biotype=true,
          pick=true 
    }

    call oncokb.oncokb as oncokbCalls {
        input:
        vcf=vep.annotatedFile,
        oncotree = oncotree,
        samplename = ctrlName,
        token = OncoKBtoken,
        searchby = "HGVSp_Short",
        AAlist = AAlist,
        canonical=canonical,
        grepRm=grepRm
    }    

    output {
        # Strelka2Germline
        Array[File]? QC_Output=normalMM_Task.picard_files
       File strelka2GermlineVCF=Strelka2Germline_Task.strelka2GermlineVCF
       # pisces outputs
       File? pisces_normal_variants=runpisces.normal_variants
       File HaplotypeVcf=CNNScoreVariantsWorkflow.cnn_filtered_vcf
       File HaplotypeVcfTbi=CNNScoreVariantsWorkflow.cnn_filtered_vcf_index
      # merged germline output
       File Merged_germline=Merge_Variants_Germline.MergedGermlineVcf
       File Merged_germlineIdx=Merge_Variants_Germline.MergedGermlineVcfIdx
       ## CNN
       File vep_annot = vep.annotatedFile
       File? vep_summary_html=vep.summary_html
       File oncokbMaf=oncokbCalls.oncokbout

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
        mkdir intervalfolder
        gatk SplitIntervals -R ${refFasta} -L ${targetIntervals} --scatter-count ${nWay} -O intervalfolder
        cp intervalfolder/*.interval_list .


    }

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        memory         : "1 GB"
    }

    output {
        Array[File] interval_files=glob("*.interval_list")
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
    String ctrlName
    
    # RUNTIME INPUT PARAMS
    String preemptible = "2"
    String diskGB_boot = "10"
    String diskGB_buffer ="5"
    String memoryGB ="4"
    String cpu ="1"
    Int minCallers = 2

    }

    # DEFAULT VALUES
    Int minV = minCallers - 1
       
    Int diskGB = ceil(size(Haplotype, "G"))*4  +  diskGB_buffer
    Boolean runS2 =if defined(STRELKA2) then true else false

    command <<<
        PISCES_pass="~{ctrlName}.Pisces.pass.vcf"
        STRELKA_pass="~{ctrlName}.S2.PASS.vcf"
        #STRELKA_unzip="~{ctrlName}.S2.unzip.vcf"
        #HP_unzip="~{ctrlName}.haplo.vcf"
        HP_pass="~{ctrlName}.haplo.pass.vcf"
        MERGED_VCF="~{ctrlName}.P_S2_HP.merged.vcf.gz"
        RENAME_MERGED_VCF_ALL="~{ctrlName}.P_S2_HP.mergedGermline.all.vcf.gz"
        RENAME_MERGED_VCF_ANN="~{ctrlName}.P_S2_HP.mergedGermline.ann.vcf"
        RENAME_MERGED_VCF_FILT="~{ctrlName}.P_S2_HP.mergedGermline.filt.vcf"


        ## filter out passed germline variants
        bcftools view -f PASS "~{STRELKA2}" > $STRELKA_pass
        bcftools view -f PASS "~{PISCES_NORMAL}" > $PISCES_pass
        bcftools view -f PASS "~{Haplotype}" > $HP_pass

        sed -i 's/##FORMAT=<ID=AD,Number=R,/##FORMAT=<ID=AD,Number=.,/g' $HP_pass

        ##bgzip $STRELKA_pass
        bgzip $PISCES_pass
        bgzip $HP_pass
        bgzip $STRELKA_pass
        tabix -p vcf $PISCES_pass.gz
        tabix -p vcf $STRELKA_pass.gz
        #bgzip $HP_pass
        tabix -p vcf $HP_pass.gz

        #merge vcfs
        bcftools merge $PISCES_pass.gz $STRELKA_pass.gz $HP_pass.gz -O vcf -o $MERGED_VCF --force-samples
        echo -e "~{ctrlName}.Pisces\n~{ctrlName}.Strelka\n~{ctrlName}.Haplotype\n" > samples.txt
        bcftools reheader -s samples.txt $MERGED_VCF > $RENAME_MERGED_VCF_ALL
        tabix -p vcf $RENAME_MERGED_VCF_ALL

        #merge
        bcftools query --format '%CHROM\t%POS\t%POS\n' $RENAME_MERGED_VCF_ALL > test.output 
        bcftools query -f '[\t%SAMPLE=%GT]\n' $RENAME_MERGED_VCF_ALL | awk '{print 3-gsub(/.\/\./, "")}' > output
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
        File MergedGermlineVcf="~{ctrlName}.P_S2_HP.mergedGermline.filt.vcf.gz"
        File MergedGermlineVcfIdx="~{ctrlName}.P_S2_HP.mergedGermline.filt.vcf.gz.tbi"
    }     
 
}

