version 1.0

task svabaCall{

    input {
    File queryBAM
    File queryIndex
    File? normalBAM
    File? normalIndex
    File reference
    File refFAIndex
    File refBWTIndex
    File refSAIndex
    File refANNIndex
    File refAMBIndex
    File refPACIndex
    String id
    Int threads = 2
    File? regions
    File? dbSNPVCF
    Int sizeAll
    Int memoryin = 32
    String runMode
    }

    String Germline = if (runMode=="Germline") then "1" else "0"

    runtime {
        docker : "quay.io/biocontainers/svaba:1.1.0--h7d7f7ad_2"
        memory :  memoryin + " GB"
        cpu : "${threads}"
        disks : "local-disk " + sizeAll + " HDD"
    }

    command {

        set -e

        if [ ${Germline} -eq "1" ]; 
        then
            germline_mode=" -L 6 -I"
        else
            germline_mode=""
        fi

        svaba run -p ${threads} \
            -G ${reference} -a ${id} ~{"-n " + normalBAM} ~{"-t " + queryBAM} ~{"-k " + regions} ~{"-D " + dbSNPVCF} --hp $germline_mode

        ls *.vcf

        gzip ~{id}.svaba.indel.vcf
        gzip ~{id}.svaba.sv.vcf
        gzip ~{id}.svaba.unfiltered.indel.vcf
        gzip ~{id}.svaba.unfiltered.sv.vcf
    }

    output {
        File? outputlog = "~{id}.log"
        File? SvABA_Indel_VCF = "~{id}.svaba.indel.vcf.gz"
        File? SvABA_SV_VCF = "~{id}.svaba.sv.vcf.gz"
        File? SvABA_Unfiltered_indel_VCF = "~{id}.svaba.unfiltered.indel.vcf.gz"
        File? SvABA_Unfiltered_SV_VCF = "~{id}.svaba.unfiltered.sv.vcf.gz"
    }

}

workflow svabaSomatic {
    input {
    File queryBAM
    File queryIndex

    File? normalBAM
    File? normalIndex

    File reference
    File refFAIndex
    File refBWTIndex
    File refSAIndex
    File refANNIndex
    File refAMBIndex
    File refPACIndex

    String id
    Int threads
    File? regions
    File? dbSNPVCF
    String runMode = "Germline" ## options: TumOnly, Paired, Germline
    }

    Int ref_size = ceil(size(reference, "GB") + size(refPACIndex, "GB") + size(refBWTIndex, "GB") + size(refSAIndex, "GB"))
    Int tumor_reads_size = ceil(size(queryBAM, "GB") + size(queryIndex, "GB"))
    Int vcf_size = if defined(dbSNPVCF) then ceil(size(dbSNPVCF, "GB")) else 0
    Int normal_reads_size = if defined(normalBAM) then ceil(size(normalBAM, "GB") + size(normalIndex, "GB")) else 0

    Int sizeAll=8*(ref_size+tumor_reads_size+vcf_size+normal_reads_size)


    call svabaCall {
        input:
            queryBAM=queryBAM,
            queryIndex=queryIndex,
            normalBAM=normalBAM,
            normalIndex=normalIndex,
            reference=reference,
            refFAIndex=refFAIndex,
            refBWTIndex=refBWTIndex,
            refSAIndex=refSAIndex,
            refANNIndex=refANNIndex,
            refAMBIndex=refAMBIndex,
            refPACIndex=refPACIndex,
            id=id,
            threads=threads,
            regions=regions,
            dbSNPVCF=dbSNPVCF,
            sizeAll=sizeAll,
            runMode=runMode
    }

    output {
        File? SvABA_Indel_VCF = svabaCall.SvABA_Indel_VCF
        File? SvABA_SV_VCF = svabaCall.SvABA_SV_VCF
        File? SvABA_Unfiltered_indel_VCF = svabaCall.SvABA_Unfiltered_indel_VCF
        File? SvABA_Unfiltered_SV_VCF = svabaCall.SvABA_Unfiltered_SV_VCF

    }

}