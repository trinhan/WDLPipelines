## This WDL pipeline runs star from single and paired-end fastq files.
## See https://github.com/alexdobin/STAR for more details
##
## 1. (Optional) Run fastqc (v0.11.8) to check quality of fastqs.
## 2. (Optional) Run Trimmomatic (v0.39) to remove low-quality basepairs and remove primers/adapter sequence.
## 3. Run STAR (v0.46.1) to align RNAseq data
##
## Required Inputs:
## - prefix: sample name
## - fastq1. Compression state doesn't matter
## - fastq2. Compression state doesn't matter
## - star_index. STAR index file, can be found in kccg-comb-bio space 
##
## Run time parameters:
## - memory: 40
## - disk_space: depends on sample, 100 should be sufficient
## - num_threads: 10
## - preemptibles: turned off, set to 0
##
## Optional Inputs: are set to default values in the STAR user guide
## - chimOutJunctionFormat: 1 Formats the output chimneric junction file
## - quantMode: set to "GeneCounts"
##
## Expected Outputs:
## - bam_file (Aligned.sortedByCoord.out.bam)
## - bam_index (Aligned.sortedByCoord.out.bam.bai)
## - chimeric_junctions (Chimeric.out.junction.gz)
## - read_counts (ReadsPerGene.out.tab.gz)
## - junctions (SJ.out.tab.gz)
## - Alignmentlogs 
## Note that chimeric bam files and transcriptome aligned bam files have been commented out

version 1.0

task star {
    input {
    File fastq1
    File? fastq2
    String prefix
    File star_index
    File? annotation_gtf

    # STAR options
    Int? outFilterMultimapNmax
    String? alignInsertionFlush
    Int? alignSplicedMateMapLminOverLmate
    Int? alignSJoverhangMin
    Int? alignSJDBoverhangMin
    Int? alignIntronMin
    Int? alignIntronMax
    Int? alignMatesGapMax
    String? alignSoftClipAtReferenceEnds
    String? alignSJstitchMismatchNmax
    Int? alignSplicedMateMapLmin
    String? chimMultimapScoreRange
    Int? chimSegmentMin
    Int? chimJunctionOverhangMin
    String? chimOutType
    Int? chimMainSegmentMultNmax
    Int? chimOutJunctionFormat
    Int? chimScoreJunctionNonGTAG
    Int? chimMultimapNmax
    Int? chimNonchimScoreDropMin
    Int? limitSjdbInsertNsj
    Int? outFilterMismatchNmax
    Float? outFilterMismatchNoverLmax
     String? outFilterType
    Float? outFilterScoreMinOverLread
    Float? outFilterMatchNminOverLread
    String? outSAMstrandField
    String? outFilterIntronMotifs
    String? outReadsUnmapped
    String? outSAMUnmapped
    String? quantMode
    String? outSAMattrRGline
    String? outSAMattributes
    Int? peOverlapNbasesMin
    Float? peOverlapMMp
    File? varVCFfile
    String? waspOutputMode
    File? sjdbFileChrStartEnd
    String readFilesCommand = "zcat"





    String? docker
    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
    Boolean runTwoPass = if defined(sjdbFileChrStartEnd) then false else true
}


    command {
        set -euo pipefail

        if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]]; then
            tar -xvvf ${fastq1}
            fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            fastq2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            if [[ $fastq1_abs == *"*_1.fastq*" ]]; then  # no paired-end FASTQs found; check for single-end FASTQ
                fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
                fastq2_abs=''
            fi
        else
            # make sure paths are absolute
            fastq1_abs=${fastq1}
            fastq2_abs=${fastq2}
            if [[ $fastq1_abs != /* ]]; then
                fastq1_abs=$PWD/$fastq1_abs
                fastq2_abs=$PWD/$fastq2_abs
            fi
        fi

        echo "FASTQs:"
        echo $fastq1_abs
        echo $fastq2_abs

        # extract index
        echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1

        #mkdir star_out
        # placeholders for optional outputs
        #touch star_out/${prefix}.Aligned.toTranscriptome.out.bam
        #touch star_out/${prefix}.Chimeric.out.sorted.bam
        #touch star_out/${prefix}.Chimeric.out.sorted.bam.bai
        #touch star_out/${prefix}.ReadsPerGene.out.tab  # run_STAR.py will gzip

        STAR \
            --genomeDir star_index \
            --readFilesIn $fastq1_abs $fastq2_abs \
            --outFileNamePrefix "${prefix}." \
            ${"--alignSJstitchMismatchNmax " + alignSJstitchMismatchNmax} \ 
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \
            ${"--alignInsertionFlush " + alignInsertionFlush} \
            ${"--alignSplicedMateMapLminOverLmate " + alignSplicedMateMapLminOverLmate} \
            ${"--alignSplicedMateMapLmin " + alignSplicedMateMapLmin} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--chimSegmentMin " + chimSegmentMin} \
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--chimOutJunctionFormat " + chimOutJunctionFormat} \
            ${"--chimMultimapScoreRange " + chimMultimapScoreRange} \
            ${"--chimScoreJunctionNonGTAG " + chimScoreJunctionNonGTAG} \
            ${"--chimMultimapNmax " + chimMultimapNmax} \
            ${"--chimNonchimScoreDropMin " + chimNonchimScoreDropMin} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            --outSAMtype BAM SortedByCoordinate \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--outSAMstrandField " + outSAMstrandField} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--outFilterType " + outFilterType} \
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--outSAMUnmapped " + outSAMUnmapped} \
            ${"--outReadsUnmapped " + outReadsUnmapped} \
            ${"--peOverlapNbasesMin " + peOverlapNbasesMin} \
            ${"--peOverlapMMp  " + peOverlapMMp} \
            ${"--quantMode " + quantMode} \
            ${"--readFilesCommand " + readFilesCommand} \
            --runThreadN ${num_threads} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            ${"--sjdbGTFfile " + annotation_gtf} \
            ${true="--twopassMode Basic" false="" runTwoPass} \
            ${"--varVCFfile " + varVCFfile} \
            ${"--waspOutputMode " + waspOutputMode} 



        ls .

        # post processing
        samtools index ${prefix}.Aligned.sortedByCoord.out.bam
        gzip ${prefix}.SJ.out.tab
        gzip ${prefix}.Chimeric.out.junction
        gzip ${prefix}.ReadsPerGene.out.tab
        ###samtools sort  --threads ${num_threads} ${prefix}.Chimeric.out.sam -o ${prefix}.Chimeric.out.sorted.bam
        ###samtools index ${prefix}.Chimeric.out.sorted.bam
    }

    output {
        File bam_file = "${prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "${prefix}.Aligned.sortedByCoord.out.bam.bai"
        ##File transcriptome_bam = "star_out/${prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "${prefix}.Chimeric.out.junction.gz"
        ##File chimeric_bam_file = "star_out/${prefix}.Chimeric.out.sorted.bam"
        ##File chimeric_bam_index = "star_out/${prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "${prefix}.ReadsPerGene.out.tab.gz"
        File junctions = "${prefix}.SJ.out.tab.gz"
        ##File junctions_pass1 = "${prefix}._STARpass1/${prefix}.SJ.pass1.out.tab.gz"
        File Finallog = "${prefix}.Log.final.out"
    }

    runtime {
        docker: select_first([docker, "broadinstitute/gtex_rnaseq:V10"])
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow star_workflow {

    call star
    
    output {
        File bam_file=star.bam_file
        File bam_index=star.bam_index
        File chimeric_junctions=star.chimeric_junctions
        File read_counts=star.read_counts
        File splice_junctions=star.junctions
        File log=star.Finallog
    }
}
