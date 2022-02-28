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
    File Bam
    File Bai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    String runMode
    File? blacklist
    File? pon
    String? preemptible
    File? bedRegion
    Int? nthreads
    Int? mem
   }

    if (runMode == "Paired" || runMode == "TumOnly"){
     call rungridssSomatic {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaFai,
            refFastaDict=refFastaDict,
            refFastaSA=refFastaSA,
            refFastaAMB=refFastaAMB,
            refFastaANN=refFastaANN,
            refFastaPAC=refFastaPAC,
            refFastaBWT=refFastaBWT,
            tumorBam=Bam,
            normalBam=normalBam,
            normalBai=normalBai,
            tumorBai=Bai,
            pairName=pairName,
            interval=interval,
            runMode=runMode,
            blacklist=blacklist,
            pondir=pon,
            preemptible=preemptible,
            bedRegion=bedRegion,
            nthreads=nthreads,
            mem=mem
        }
    }


    if (runMode == "Germline"){
     call rungridssGermline {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaFai,
            refFastaDict=refFastaDict,
            refFastaSA=refFastaSA,
            refFastaAMB=refFastaAMB,
            refFastaANN=refFastaANN,
            refFastaPAC=refFastaPAC,
            refFastaBWT=refFastaBWT,
            bam=Bam,
            bai=Bai,
            pairName=pairName,
            interval=interval,
            runMode=runMode,
            blacklist=blacklist,
            preemptible=preemptible,
            bedRegion=bedRegion,
            nthreads=nthreads,
            mem=mem
        }
    }

    output {
        File? somatic_high_conf = rungridssSomatic.somatic_high_conf
        File? somatic_all = rungridssSomatic.somatic_all
        File? germlineOutputs = rungridssGermline.Germline
        File? supportingBamSomatic = rungridssSomatic.evidence_bam
        File? supportingBamGerm = rungridssGermline.evidence_bam

    }
}

task rungridssSomatic {
    input {
    File refFasta
    File refFastaFai
    File refFastaDict
    File refFastaSA
    File refFastaAMB
    File refFastaANN
    File refFastaPAC
    File refFastaBWT
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    String runMode
    File? blacklist
    File? pondir
    File? bedRegion
    }

    String normP = select_first([ normalBam ,""])
    String tumP = select_first([ tumorBam ,""])
    String tumPrefix=if (tumP!="") then basename(sub(tumP,"\\.bam$", "")) else ""
    String normPrefix= if (normP!="") then basename(sub(normP,"\\.bam$", "")) else ""

    String matchPair = if (runMode=="Paired" ||  ( defined(normalBam) && defined(tumorBam))) then "1" else "0"
    Int tumBamSize = select_first([ceil(size(tumorBam, "G")+size(tumorBai, "G")), 0])
    Int normBamSize = select_first([ceil(size(normalBam, "G")+size(normalBai, "G")), 0])
    Int disk_size=2*ceil(tumBamSize+ normBamSize+3)
   
    command <<<
        set -e

        ls /opt/gridss
        mkdir pondir 
        tar -C pondir  -xvf ~{pondir} 

        GRIDSS_JAR=/opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar
        ##PATH=/opt/gridss/:/opt/RepeatMasker:/opt/rmblast/:/opt/trf:/opt/kraken2:/opt/blast:/opt/edirect:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
        echo $GRIDSS_JAR
        # define how to run normal if present
        runCommand=""
        if [ ~{matchPair} -eq "1" ]; then
        runCommand=" ~{normalBam} "
        fi

        # run gridss
        gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -o ~{pairName}.vcf ~{"--targetbed " + bedRegion} \
        ~{"--blacklist "+ blacklist} --threads ~{nthreads} $runCommand ~{tumorBam}
  
        # run somatic filter
        gridss_somatic_filter --pondir pondir --input ~{pairName}.vcf --output ~{pairName}.high_confidence_somatic.vcf.gz \
        --fulloutput ~{pairName}.high_and_low_confidence_somatic.vcf.gz --scriptdir $(dirname $(which gridss_somatic_filter)) -n 1 -t 2

        ls *.bam*

    >>>

    output {
        File somatic_high_conf= "~{pairName}.high_confidence_somatic.vcf.gz" 
        File somatic_all= "~{pairName}.high_and_low_confidence_somatic.vcf.gz" 
        File evidence_bam= "~{pairName}.vcf.assembly.bam"

    }

    runtime {
        docker: "gridss/gridss:2.13.1"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}

task rungridssGermline {
    input {
    File refFasta
    File refFastaFai
    File refFastaDict
    File refFastaSA
    File refFastaAMB
    File refFastaANN
    File refFastaPAC
    File refFastaBWT
    File bam
    File bai
    String pairName
    File? interval
    Int? nthreads = 3
    String mem = 14
    Int preemptible =3
    String runMode
    File? blacklist
    File? bedRegion
    
    }

    Int disk_size=2*ceil((size(bam , "GB")+size(refFasta, "GB")+3))
   
    command <<<
        set -e

        ls /opt/gridss
        ####mkdir pondir 
        ######tar -C pondir -xvf ${pondir} 

        GRIDSS_JAR=/opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar
        ##PATH=/opt/gridss/:/opt/RepeatMasker:/opt/rmblast/:/opt/trf:/opt/kraken2:/opt/blast:/opt/edirect:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
        echo $GRIDSS_JAR
        # define how to run normal if present

        gridss -r ~{refFasta} --jar /opt/gridss/gridss-2.13.1-gridss-jar-with-dependencies.jar -o ~{pairName}.Germline.vcf  ~{"--targetbed " + bedRegion} \
        ~{"--blacklist " + blacklist} --threads ~{nthreads} ~{bam}

        mv ~{pairName}.Germline.vcf.assembly.bam.gridss.working/~{pairName}.Germline.vcf.assembly.bam.sv.bam .
        mv ~{pairName}.Germline.vcf.assembly.bam.gridss.working/~{pairName}.Germline.vcf.assembly.bam.sv.bam.bai .

        ls *

    >>>

    output {
        File Germline= "~{pairName}.Germline.vcf" 
        File evidence_bam= "~{pairName}.Germline.vcf.assembly.bam.sv.bam"
        File evidence_bai= "~{pairName}.Germline.vcf.assembly.bam.sv.bam.bai"
        Array[File] log= glob("*.log") 

    }

    runtime {
        docker: "gridss/gridss:2.13.1"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}
