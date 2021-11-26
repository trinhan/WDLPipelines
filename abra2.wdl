## ABRA2 realigner

version 1.0

workflow abra2 {
    input {
 File vcf
 File? vcfidx
 String sample_name
 File? normalBam
 File? normalIdx
File tumorBam
File tumorIdx
File refFasta
File refFastaIdx
File? targets
} 


call runabra2 {
        input:
        normalBam=normalBam,
        normalIdx=normalIdx,
        tumorBam=tumorBam,
        tumorIdx=tumorIdx,
        refFasta=refFasta,
        vcf=vcf,
        vcfidx=vcfidx,
        sample_name=sample_name,
        refFastaIdx=refFastaIdx,
        targets=targets
    }
}

task runabra2 {
    input {
    File? normalBam
    File? normalIdx
    String sample_name
    File tumorBam
    File tumorIdx
    File refFasta
    File refFastaIdx
    File vcf
    File? vcfidx
    File? targets
    
    Int cpu = 3
    Int machine_mem_gb = 12
    Int? preemptible_tries =3
    Int? max_retries
     }

        Int command_mem=machine_mem_gb-1
        Int disk_space_gb=4*ceil(size(tumorBam, "GB")+size(normalBam, "GB")+size(refFasta, "GB"))
        String existBam= if defined(normalBam) then "1" else "0"


    command {

        set -e

        mkdir tmp

        if [ ${existBam} -eq "1" ];
        then
        echo 'normal exists'
                java -Xmx${command_mem}g -jar /usr/local/bin/abra2.jar --in ${normalBam},${tumorBam} --out ${sample_name}.N.abra2.bam,${sample_name}.T.abra2.bam \
                --ref ${refFasta} --threads ${cpu} --in-vcf ${vcf} ${"--targets " + targets} --tmpdir tmp > ${sample_name}.log

                cp "${sample_name}.N.abra2.bai" "${sample_name}.N.abra2.bam.bai" 
        else
        echo 'no normal'
                java -Xmx${command_mem}g -jar /usr/local/bin/abra2.jar --in ${tumorBam} --out ${sample_name}.T.abra2.bam \
                --ref ${refFasta} --threads ${cpu} --in-vcf ${vcf} ${"--targets " + targets} --tmpdir tmp > ${sample_name}.log

                cp "${sample_name}.T.abra2.bai" "${sample_name}.T.abra2.bam.bai" 
        fi
       
    }

    runtime {
          docker: "mskaccess/abra2:2.22"
          preemptible: select_first([preemptible_tries, 1])
          maxRetries: select_first([max_retries, 1])
          memory: machine_mem_gb + " GB"
          disks: "local-disk " + disk_space_gb + " HDD"
          cpu: cpu
    }
    output {
        File Abra2Tumour = "${sample_name}.T.abra2.bam"
        File Abra2TumourIdx = "${sample_name}.T.abra2.bam.bai"
        File? Abra2Normal = "${sample_name}.N.abra2.bam"
        File? Abra2NormalIdx = "${sample_name}.N.abra2.bam.bai"
        }

}
