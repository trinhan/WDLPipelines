## This WDL pipeline runs RSeQC (version4.0.0) from a bam file
## See http://rseqc.sourceforge.net/#download-rseqc for more information
##
## Run options: The following modules are supported
## run_read_dup: read_duplication.py   
## run_gene_body: geneBody_coverage.py
## run_fpkm_uq: FPKM-UQ.py 
## run_bam_stat: bam_stat.py 
## run_read_dist: read_distribution.py
##
## Required Inputs:
## - sampleName: sample name
## - bam: sample bam
## - bam_idx: sample bam idx
##
## Annotation Files:
## These can be accessed here: http://rseqc.sourceforge.net/#download-rseqc
##    File? refBed - reference transcriptome bed for read_distribution 
##    File? housekeepBed - house keeping genes for gene body coverage analysis
##    File? gencodeannotationGtf - see http://rseqc.sourceforge.net/#fpkm-uq-py: for FPKM similar to TCGA
##    File? gencodegeneinfoTsv - see http://rseqc.sourceforge.net/#fpkm-uq-py: for FPKM similar to TCGA
## 
## Expected Outputs:
## - read_duplicates 
## - genebody 
## - FPKM_UQ
## - run_read_dist 
## - run_bam_stat

version 1.0

workflow RSeQC {
    call RunQCChecks
    output {
        Array[File]? read_duplicates = RunQCChecks.read_duplicates
        Array[File]? genebody = RunQCChecks.genebody
        Array[File]? FPKM_UQ = RunQCChecks.FPKM_UQ
        File? run_read_dist = RunQCChecks.run_read_dist
        File? run_bam_stat = RunQCChecks.run_bam_stat
    }
}

task RunQCChecks {
    input {
    File bam
    File bamidx
    String sampleName
    String run_read_dup
    String run_gene_body
    String run_fpkm_uq
    String run_bam_stat
    String run_read_dist
    File? refBed
    File? housekeepBed
    File? gencodeannotationGtf
    File? gencodegeneinfoTsv
}

    Int diskspace = ceil(3*size(bam, "GB"))

command {

    if [ "${run_read_dup}" = true ] ;
        then 
        echo "Run read duplication"
        read_duplication.py -i ~{bam} -o ~{sampleName}_read_duplicates
    fi

    if [ "${run_gene_body}" = true ] ;
        then 
        echo "Run gene body"
        geneBody_coverage.py -r ~{housekeepBed} -i ~{bam}  -o ~{sampleName}_genebody
    fi

    if [ "${run_fpkm_uq}" = true ] ;
        then
        echo "Run fpkm"
        FPKM-UQ.py --bam ~{bam} --gtf ~{gencodeannotationGtf} --info ~{gencodegeneinfoTsv} -o ~{sampleName}_fpkm
    fi

    if [ "${run_bam_stat}" = true ] ;
        then 
        echo "Run bam stat"
        bam_stat.py -i ~{bam} > ~{sampleName}.bam_stat
    fi 

    if [ "${run_read_dist}" = true ] ;
        then 
        echo "Run read read_distribution"
        if [[ ${refBed} == *".gz" ]];
            then 
            gunzip -c ~{refBed} > refBed.bed
            read_distribution.py  -i ~{bam} -r refBed.bed > ~{sampleName}.read_distribution
            else
            read_distribution.py  -i ~{bam} -r ~{refBed} > ~{sampleName}.read_distribution
        fi 
    fi
}

output {
    Array[File]? read_duplicates = glob("~{sampleName}_read_duplicates*")
    Array[File]? genebody = glob("~{sampleName}_genebody*")
    Array[File]? FPKM_UQ = "~{sampleName}_fpkm*"
    File? run_read_dist = "~{sampleName}.read_distribution"
    File? run_bam_stat = "~{sampleName}.bam_stat"
}

runtime { 
    docker: "trinhanne/rseqc:4.0.0"
    preemptible: "2"
    memory: "5 GB"
    disks: "local-disk ~{diskspace} HDD"
}

}