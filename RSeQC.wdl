## This WDL pipeline runs RSeQC (version4.0.0) from a bam file
## See http://rseqc.sourceforge.net/#download-rseqc for more information
##
## 1. (Optional) Run fastqc (v0.11.8) to check quality of fastqs.
## 2. (Optional) Run Trimmomatic (v0.39) to remove low-quality basepairs and remove primers/adapter sequence.
## 3. Run STAR (v0.46.1) to align RNAseq data
##
## Required Inputs:
## - prefix: sample name
## - bam: sample bam
## - bam_idx: sample bam idx
## - star_index. STAR index file, can be found in kccg-comb-bio space 
##
## Modules set up:
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

workflow RSeQC {
    call RunQCChecks
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
        read_duplication.py -i ~{bam} -o ~{sampleName}_read_duplicates
    fi

    if [ "${run_gene_body}" = true ] ;
        then 
    geneBody_coverage.py -r ~{housekeepBed} -i ~{bam}  -o ~{sampleName}_genebody
    fi

    if [ "${run_fpkm_uq}" = true ] ;
        then
    FPKM-UQ.py --bam ~{bam} --gtf ~{gencodeannotationGtf} --info ~{gencodegeneinfoTsv} -o ~{sampleName}_fpkm
    fi

    if [ "${run_bam_stat}" = true ] ;
        then 
        bam_stat.py -i ~{bam} > ~{sampleName}.bam_stat
    fi 

    if [ "${run_read_dist}" = true ] ;
        then 
        read_distribution.py  -i ~{bam} -r ~{refBed} > ~{sampleName}.read_distribution
    fi
}

output {
    Array[File] read_duplicates = glob("~{sampleName}_read_duplicates")
    Array[File] genebody = glob("~{sampleName}_genebody")
    File FPKM_UQ = "~{sampleName}_fpkm"
    File run_read_dist = "~{sampleName}.read_distribution"
    File run_bam_stat = "~{sampleName}.bam_stat"
}

runtime { 
    docker: "clinicalgenomics/rseqc:4.0.0"
    preemptible: "2"
    memory: "5 GB"
    disks: "local-disk ~{diskspace} HDD"
}

}