## Run blat

version 1.0


task blat {
    input{
    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File MAF
    File genome_bit
    String pairName

    # FILE SIZE
    Int tumorBam_size

    # RUNTIME INPUT PARAMS
    String preemptible ="1"
    String diskGB_boot = "15"
    String diskGB_buffer = "15"
    String machine_memoryGB = "10"
    String cpu ="1"
    ## can also be saved as /opt/hg38.2bit
    String savPath = "/opt/hg19.2bit" 
}

    # COMPUTE DISK SIZE
    Int diskGB = 2*ceil(tumorBam_size + size(MAF, "G") + size(genome_bit, "G") + diskGB_buffer)

    parameter_meta {
        tumorBam : "sample tumor BAM file"
        tumorBamIdx : "sample tumor BAI file (indexed BAM)"
        MAF : "filename pointing to a mutation annotation format (MAF) file (data for somatic point mutations)"
        pairName : "tumor sample name, string prefix of the output"
    }

    command {

        set -euxo pipefail

        cp -v ${genome_bit} ${savPath}
        python /opt/realign.py ${tumorBam} ${MAF} ${pairName}

        # Count number of passed and rejected mutations
        python /usr/local/bin/count_variants.py "${pairName}.blat.maf" "${pairName}.count_passed_mutations.txt"
        python /usr/local/bin/count_variants.py "${pairName}.blat.rejected.maf" "${pairName}.count_rejected_mutations.txt"

        python /usr/local/bin/add_judgement_column.py \
        --input "${pairName}.blat.all.maf" \
        --output "${pairName}.blat.all.with_judgement_annotations.maf" \
        --column "realign_judgment" \
        --pass_flag "KEEP"

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
        Int num_rejected_mutations=read_int("${pairName}.count_rejected_mutations.txt")
        Int num_passed_mutations=read_int("${pairName}.count_passed_mutations.txt")
        File passMaf="${pairName}.blat.maf"
        File debug_results="${pairName}.blat.rejected.maf"
        File allMaf="${pairName}.blat.all.with_judgement_annotations.maf"
    }
}

workflow runBlat{
   input{
    # TASK INPUT PARAMS
    File tumorBam
    File tumorBamIdx
    File MAF
    File genome_bit
    String pairName

    # FILE SIZE
    Int tumorBam_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_boot
    String diskGB_buffer
    String memoryGB
    String cpu
}

	call blat {
        input:
            tumorBam=tumorBam,
            tumorBamIdx=tumorBamIdx,
            MAF=MAF, 
            pairName=pairName,
            genome_bit=genome_bit,
            tumorBam_size=tumorBam_size
        }

}
