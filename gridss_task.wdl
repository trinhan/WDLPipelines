version 1.0

workflow gridss_workflow {
  input {
    File refFasta
    File refFastaFai
    File refFastaDict
    File refFastaSA
    File refFastaAMB
    File refFastaANN
    File refFastaPAC
    File refFastaBWT
    File? tumorBam
    File? tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    String runMode
    File? blacklist
   }

     call rungridss {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaFai,
            refFastaDict=refFastaDict,
            refFastaSA=refFastaSA,
            refFastaAMB=refFastaAMB,
            refFastaANN=refFastaANN,
            refFastaPAC=refFastaPAC,
            refFastaBWT=refFastaBWT,
            tumorBam=tumorBam,
            normalBam=normalBam,
            normalBai=normalBai,
            tumorBai=tumorBai,
            pairName=pairName,
            interval=interval,
            runMode=runMode
        }

    output {
        File? tumor_unique_variants_phased=rungridss.tumor_unique_variants_phased
        File? tumor_unique_variants=rungridss.tumor_unique_variants
        File? normal_variants_same_site=rungridss.normal_variants_same_site
        File? normal_variants=rungridss.normal_variants
        File? tumor_variants=rungridss.tumor_variants
        File? venn_zip=rungridss.venn_zip

    }
}

task rungridss {
    input {
    File refFasta
    File refFastaFai
    File refFastaDict
    File refFastaSA
    File refFastaAMB
    File refFastaANN
    File refFastaPAC
    File refFastaBWT
    File? tumorBam
    File? tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    String runMode
    }


    String normP = select_first([ normalBam ,""])
    String tumP = select_first([ tumorBam ,""])
    String tumPrefix=if (tumP!="") then basename(sub(tumP,"\\.bam$", "")) else ""
    String normPrefix= if (normP!="") then basename(sub(normP,"\\.bam$", "")) else ""

    String runTum = if (runMode=="TumOnly") then "1" else "0"
    String runGerm = if (runMode=="Germline" ||  defined(normalBam)) then "1" else "0"
    String matchPair = if (runMode=="Paired" ||  ( defined(normalBam) && defined(tumorBam))) then "1" else "0"
    Int tumBamSize = select_first([ceil(size(tumorBam, "G")+size(tumorBai, "G")), 0])
    Int normBamSize = select_first([ceil(size(normalBam, "G")+size(normalBai, "G")), 0])
    Int disk_size=2*(tumBamSize+ normBamSize+3)
   
    command <<<
        set -e

        ls /opt/gridss

        ##GRIDSS_JAR=/opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar
        ##PATH=/opt/gridss/:/opt/RepeatMasker:/opt/rmblast/:/opt/trf:/opt/kraken2:/opt/blast:/opt/edirect:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

        echo $GRIDSS_JAR

        if [ ~{runTum} -eq "1" ];
        then        
        ## run preprocessing
        gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -s preprocess -t 4 ~{tumorBam}
        fi
        
        ## Run germline based calling

        if [ ~{runGerm} -eq "1" ];
        then
        gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -s preprocess -t 4 ~{normalBam}
        ## dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf somatic_~{pairName}/~{normPrefix}.recal.vcf --bam ~{normalBam}
        fi
        
        if [ ~{matchPair} -eq "1" ];
        then 
            gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -s preprocess -t 4 ~{tumorBam}
            gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -s preprocess -t 4 ~{normalBam}      
            gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -t 8 -s assemble -a assembly.bam ~{normalBam} ~{tumorBam}
        fi

        
    >>>

    output {
        File? tumor_unique_variants= "~{tumPrefix}.somatic.unique.recal.vcf" 
        File? tumor_unique_variants_phased= "~{tumPrefix}.somatic.unique.recal.phased.vcf" 
        File? normal_variants_same_site= "variant2_~{pairName}/~{normPrefix}.genome.recal.vcf" 
        File? normal_variants_recal =  "somatic_${pairName}/~{normPrefix}.recal.vcf" 
        File? normal_variants =  "somatic_~{pairName}/~{normPrefix}.vcf" 
        File? tumor_variants =  "somatic_${pairName}/~{tumPrefix}.recal.vcf" 
        File? venn_zip="~{pairName}_venn_gridss.tar.gz"


    }

    runtime {
        docker: "gridss/gridss:2.13.2"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}
