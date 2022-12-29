## Steps in this workflow:
## 1. Parallelise the intervals
## 2. run vardict on each
## 3. update headers
## 4. consolidate into 1 final vcf file

version 1.0

workflow VardictWF {
    input {
        File refFastaIdx
        File refFasta
        File refFastaDict
        File? normalBam
        File? normalBamIdx
        File tumorBam
        File tumorBamIdx
        String? ctrlName
        String caseName
        Array[File] bed_list
        Array[Int] scatterIndices 
        String gatk_docker
    }

    scatter (idx in scatterIndices) {
        call runVardict {
            input:
                referenceFasta=refFasta,
                referenceFastaFai=refFastaIdx,
                tumorBam=tumorBam,
                tumorBamIndex=tumorBamIdx,
                normalBam=normalBam,
                normalBamIndex=normalBamIdx,
                outputName=caseName,
                tumorSampleName=caseName,
                bedFile=bed_list[idx],
                tumorSampleName=caseName,
                normalSampleName=ctrlName
        }
    }

     call UpdateHeaders {
        input:
            input_vcfs = runVardict.vcfFile,
            ref_dict=refFastaDict,
            gatk_docker = gatk_docker,
            caller = "Vardict"
    }

    call CombineVariants {
        input:
            input_header=UpdateHeaders.head_vcf,
            ref_fasta = refFasta,
            ref_fai = refFastaIdx,
            ref_dict = refFastaDict,
            gatk_docker = gatk_docker,
            sample_name = caseName,
            caller="Vardict"
    }

    output {
        File vardict=CombineVariants.merged_vcf
    }

}


task runVardict {
    input {
        ## sample information here
        String tumorSampleName
        File tumorBam
        File tumorBamIndex
        String? normalSampleName
        File? normalBam
        File? normalBamIndex
        ## reference files go here
        File referenceFasta
        File referenceFastaFai
        File? bedFile
        String outputName
        Boolean outputCandidateSomaticOnly = true
        Boolean outputAllVariantsAtSamePosition = true
        Boolean callSVs = false
        Float mappingQuality = 20
        Int minimumTotalDepth = 8
        Int minimumVariantDepth = 5
        Float minimumAlleleFrequency = 0.05
        Int chromosomeColumn = 1
        Int startColumn = 2
        Int endColumn = 3
        Int geneColumn = 4
        Int runLocalRelignment = 0
        ## run time parameters
        String javaXmx = "15G"
        Int threads = 4
        String memory = "16"
        String dockerImage = "quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0"
        Int? opt_preempt
        Int? timeMinutes = 350
    }
        Int diskGB=3*ceil(size(tumorBam, "GB")+size(normalBam, "GB")+size(referenceFasta, "GB"))
        String runPaired = if defined(normalBam) then "1" else "0"

    command {


        ## remove this command:  -XX:ParallelGCThreads=1
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{javaXmx}"
        vardict-java \
        ~{"-th " + threads} \
        -G ~{referenceFasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam}~{"|" + normalBam}" \
        ~{false="--nosv " true = "" callSVs} \
        ~{true="" false="-z" defined(normalBam)} \
        -c ~{chromosomeColumn} \
        -S ~{startColumn} \
        -E ~{endColumn} \
        -g ~{geneColumn} \
        ~{bedFile} | \
        ~{true="testsomatic.R" false="teststrandbias.R" defined(normalBam)} | \
        ~{true="var2vcf_paired.pl" false="var2vcf_valid.pl" defined(normalBam)} \
        -N "~{tumorSampleName}~{"|" + normalSampleName}" \
        ~{true="" false="-E" defined(normalBam)} \
        ~{true="-M" false="" outputCandidateSomaticOnly} \
        ~{true="-A" false="" outputAllVariantsAtSamePosition} \
        -Q ~{mappingQuality} \
        -d ~{minimumTotalDepth} \
        -v ~{minimumVariantDepth} \
        -f ~{minimumAlleleFrequency} \
        > ~{outputName}.vardict.vcf
        
    }

    output {
        File vcfFile = "~{outputName}.vardict.vcf"
    }

    runtime {
        cpu: threads
        memory: "${memory} GB"
        docker: dockerImage
        disks: "local-disk ${diskGB} HDD"
        preemptible: select_first([opt_preempt, 1])
        time_minutes: timeMinutes
    }

    parameter_meta {
        # inputs
        tumorSampleName: {description: "The name of the tumor/case sample.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        bedFile: {description: "A bed file describing the regions to operate on. These regions must be below 1e6 bases in size.", category: "required"}
        outputCandidateSomaticOnly: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-M` flag.", category: "advanced"}
        outputAllVariantsAtSamePosition: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-A` flag.", category: "advanced"}
        mappingQuality: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-Q` option.", category: "advanced"}
        minimumTotalDepth: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-d` option.", category: "advanced"}
        minimumVariantDepth: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-v` option.", category: "advanced"}
        minimumAlleleFrequency: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-f` option.", category: "advanced"}
        chromosomeColumn: {description: "Equivalent to vardict-java's `-c` option.", category: "advanced"}
        startColumn: {description: "Equivalent to vardict-java's `-S` option.", category: "advanced"}
        endColumn: {description: "Equivalent to vardict-java's `-E` option.", category: "advanced"}
        geneColumn: {description: "Equivalent to vardict-java's `-g` option.", category: "advanced"}
        normalSampleName: {description: "The name of the normal/control sample.", category: "common"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "common"}
        normalBamIndex: {description: "The normal/control sample's BAM file.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}   
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