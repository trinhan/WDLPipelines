version 1.0


task VarDict {
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
        Int runLocalRelignment = 1
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

     if [ ~{runPaired} -eq "1" ]; then
        ## remove this command:  -XX:ParallelGCThreads=1
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{javaXmx}"
        vardict-java \
        ~{"-th " + threads} \
        -G ~{referenceFasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam}~{"|" + normalBam}" \
        ~{true="" false="-z" defined(normalBam)} \
        ~{false="-U" true="" callSVs} \ 
        -k ~{runLocalRelignment} \
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
    else 
        ## remove this command:  -XX:ParallelGCThreads=1
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{javaXmx}"
        vardict-java \
        ~{"-th " + threads} \
        -G ~{referenceFasta} \
        -N ~{tumorSampleName} \
        -b ~{tumorBam} \
        ~{true="" false="-z" defined(normalBam)} \
        ~{false="-U" true="" callSVs} \ 
        -k ~{runLocalRelignment} \
        -c ~{chromosomeColumn} \
        -S ~{startColumn} \
        -E ~{endColumn} \
        -g ~{geneColumn} \
        ~{bedFile} | \
        teststrandbias.R | var2vcf_valid.pl \
        -N ~{tumorSampleName} -E \
        ~{true="-M" false="" outputCandidateSomaticOnly} \
        ~{true="-A" false="" outputAllVariantsAtSamePosition} \
        -Q ~{mappingQuality} \
        -d ~{minimumTotalDepth} \
        -v ~{minimumVariantDepth} \
        -f ~{minimumAlleleFrequency} \
        > ~{outputName}.vardict.vcf
    fi 
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