version 1.0

task svabaCall{

    input {
    File tumorBAM
    File tumorIndex
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
    
    }

    runtime {
        docker : "quay.io/biocontainers/svaba:1.1.0--h7d7f7ad_2"
        memory :  memoryin + " GB"
        cpu : "${threads}"
        disks : "local-disk " + sizeAll + " HDD"
    }

    command {

        set -e

        if [[ -e ${tumorBAM} ]]; 
        then
            germline_mode=""
        else
            germline_mode=" -L 6 -I"
        fi
        
        echo $germline_mode

        svaba run -p ${threads} \
            -G ${reference} -a ${id} ~{"-n " + normalBAM} ~{"-t " + tumorBAM} ~{"-k " + regions}  ~{"-D " + dbSNPVCF} --hp $germline_mode

        ls 
          
    }

    output {
        File? outputlog = "${id}.log"
        File? SvABA_Somatic_Indel_VCF = "${id}.svaba.somatic.indel.vcf"
        File? SvABA_Somatic_SV_VCF = "${id}.svaba.somatic.sv.vcf"
        File? SvABA_Somatic_Unfiltered_indel_VCF = "${id}.svaba.unfiltered.somatic.indel.vcf"
        File? SvABA_Somatic_Unfiltered_SV_VCF = "${id}.svaba.unfiltered.somatic.sv.vcf"
        File? SvABA_Germline_Indel_VCF = "${id}.svaba.germline.indel.vcf"
        File? SvABA_Germline_SV_VCF = "${id}.svaba.germline.sv.vcf"
    }

}

workflow svabaSomatic {
    input {
    File? tumorBAM
    File? tumorIndex

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
    }

    Int ref_size = ceil(size(reference, "GB") + size(refPACIndex, "GB") + size(refBWTIndex, "GB") + size(refSAIndex, "GB"))
    Int tumor_reads_size = ceil(size(tumorBAM, "GB") + size(tumorIndex, "GB"))
    Int vcf_size = if defined(dbSNPVCF) then ceil(size(dbSNPVCF, "GB")) else 0
    Int normal_reads_size = if defined(normalBAM) then ceil(size(normalBAM, "GB") + size(normalIndex, "GB")) else 0

    Int sizeAll=8*(ref_size+tumor_reads_size+vcf_size+normal_reads_size)


    call svabaCall {
        input:
            tumorBAM=tumorBAM,
            tumorIndex=tumorIndex,
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
            sizeAll=sizeAll
    }

}