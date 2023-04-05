version 1.0

workflow pisces_workflow {
  input {
    File? refFasta
    File? refFastaFai
    File? refFastaDict
    File? pisces_reference
    File? tumorBam
    File? tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    String runMode
   }

 if (runMode =="Paired"){
     call runpiscesSomaticPaired {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaFai,
            refFastaDict=refFastaDict,
            tumorBam=tumorBam,
            normalBam=normalBam,
            normalBai=normalBai,
            tumorBai=tumorBai,
            pairName=pairName,
            pisces_reference=pisces_reference,
            interval=interval,
            runMode=runMode
        }
    }

if (runMode == "Germline"){
     call runpiscesGermline {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaFai,
            refFastaDict=refFastaDict,
            normalBam=normalBam,
            normalBai=normalBai,
            pairName=pairName,
            pisces_reference=pisces_reference,
            interval=interval,
            runMode=runMode
        }
    }

    output {
        File tumor_unique_variants_phased=select_first([runpiscesSomatic.tumor_unique_variants_phased, "NULL"])
        File tumor_unique_variants=select_first([runpiscesSomatic.tumor_unique_variants, "NULL"])
        File normal_variants=select_first([runpiscesSomatic.normal_variants_same_site, runpiscesGermline.normal_variants])
        File? venn_zip=select_first([runpiscesSomatic.venn_zip, runpiscesGermline.venn_zip])
        File? refzip=select_first([runpiscesSomatic.refzip, runpiscesGermline.refzip])

    }
}

task runpiscesGermline {
    input {
    File? refFasta
    File? refFastaFai
    File? refFastaDict
    File? pisces_reference
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    Boolean MNVcall = true
    Boolean scatterchr = true
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    String saveDict = "0"
    String runMode
    String runScylla = "0"
    }


    String normP = select_first([ normalBam ,""])
    String normPrefix= if (normP!="") then basename(sub(normP,"\\.bam$", "")) else ""
    String buildRef = if defined(pisces_reference) then "0" else "1"

    Int disk_size=2*(ceil(size(normalBam, "G")+size(normalBai, "G")+size(refFasta, "G")+size(pisces_reference, "G")))
   
    command <<<
        set -e

        mkdir somatic_~{pairName}
        sname=~{pisces_reference}

        ## Build the reference libraries here if they do not exist

        if [ ~{buildRef} -eq "1" ];
        then
        echo 'create the reference'
        sname=`basename ~{refFasta}`
        sname=$(echo $sname| cut -f 1 -d '.')
        echo $sname
        mkdir $sname
        cp ~{refFasta} $sname
        cp ~{refFastaFai} $sname
        cp ~{refFastaDict} $sname
        dotnet /app/CreateGenomeSizeFile_5.2.10.49/CreateGenomeSizeFile.dll -g $sname -s "Homo sapien $sname" -o $sname
        fi 
        
        if [ ~{buildRef} -eq "0" ];
        then
        ## to unpack the tar file
        tar xvzf ~{pisces_reference}
        ## use this when using a reference file
        ## cp -r ~{pisces_reference} . 
        sname=`basename ~{pisces_reference}`
        sname=$(echo $sname| cut -f 1 -d '.')
        fi

        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{normalBam} -CallMNVs ~{MNVcall} -gVCF false --threadbychr ~{scatterchr} --collapse true -ploidy diploid \
        --filterduplicates true --maxthreads ~{nthreads} -o somatic_~{pairName} ~{"-i " + interval}

        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_~{pairName}/~{normPrefix}.vcf --out somatic_~{pairName}

        if [[ -f somatic_~{pairName}/~{normPrefix}.vcf.recal ]];
            then
              cp somatic_~{pairName}/~{normPrefix}.vcf.recal somatic_~{pairName}/~{normPrefix}.recal.vcf
            elif [[ -f somatic_~{pairName}/~{normPrefix}.vcf ]];
            then
              cp somatic_~{pairName}/~{normPrefix}.vcf somatic_~{pairName}/~{normPrefix}.recal.vcf
        fi
        ## dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf somatic_~{pairName}/~{normPrefix}.recal.vcf --bam ~{normalBam}
        

        if [ ~{runScylla} -eq "1" ];
        then  
        ## perform phasing here
        dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf ~{normPrefix}.somatic.unique.recal.vcf --bam ~{normalBam}
        fi

        if [ ~{saveDict} -eq "1" ];
        then
            tar -zvcf refPisces.tar.gz $sname
        fi
        
    >>>

    output {
        File normal_variants =  select_first(["somatic_${pairName}/~{normPrefix}.recal.vcf", "somatic_~{pairName}/~{normPrefix}.vcf"])
        File? venn_zip="~{pairName}_venn_pisces.tar.gz"
        File? refzip="refPisces.tar.gz"

    }

    runtime {
        docker: "trinhanne/pisces:5.2.10"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}


task runpiscesSomaticPaired {
    input {
    File? refFasta
    File? refFastaFai
    File? refFastaDict
    File? pisces_reference
    File? tumorBam
    File? tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    Boolean MNVcall = true
    Boolean scatterchr = true
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    String saveDict = "0"
    String runMode
    String runScylla = "0"
    }


    String normP = select_first([ normalBam ,""])
    String tumP = select_first([ tumorBam ,""])
    String tumPrefix=if (tumP!="") then basename(sub(tumP,"\\.bam$", "")) else ""
    String normPrefix= if (normP!="") then basename(sub(normP,"\\.bam$", "")) else ""
    String buildRef = if defined(pisces_reference) then "0" else "1"

    String runTum = if (runMode!="Germline") then "1" else "0"
    String runGerm = if (runMode!="TumOnly" ||  defined(normalBam)) then "1" else "0"
    String matchPair = if (runMode=="Paired" ||  ( defined(normalBam) && defined(tumorBam))) then "1" else "0"
    String tOnly = if (runMode=="TumOnly" ||  ( !defined(normalBam) && defined(tumorBam))) then "1" else "0"
    Int tumBamSize = select_first([ceil(size(tumorBam, "G")+size(tumorBai, "G")), 0])
    Int normBamSize = select_first([ceil(size(normalBam, "G")+size(normalBai, "G")), 0])
    Int disk_size=2*(tumBamSize+ normBamSize+3)
   
    command <<<
        set -e

        mkdir somatic_~{pairName}
        sname=~{pisces_reference}

        ## Build the reference libraries here if they do not exist

        if [ ~{buildRef} -eq "1" ];
        then
        echo 'create the reference'
        sname=`basename ~{refFasta}`
        sname=$(echo $sname| cut -f 1 -d '.')
        echo $sname
        mkdir $sname
        cp ~{refFasta} $sname
        cp ~{refFastaFai} $sname
        cp ~{refFastaDict} $sname
        dotnet /app/CreateGenomeSizeFile_5.2.10.49/CreateGenomeSizeFile.dll -g $sname -s "Homo sapien $sname" -o $sname
        fi 
        
        if [ ~{buildRef} -eq "0" ];
        then
        ## to unpack the tar file
        tar xvzf ~{pisces_reference}
        ## use this when using a reference file
        ## cp -r ~{pisces_reference} . 
        sname=`basename ~{pisces_reference}`
        sname=$(echo $sname| cut -f 1 -d '.')
        fi

        ## Run the tumor based calling of variants

        if [ ~{runTum} -eq "1" ];
        then        
        ## run the variant calling: this works
        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{tumorBam} -CallMNVs ~{MNVcall} -gVCF false --threadbychr ~{scatterchr} --collapse true -c 5 \
        --filterduplicates true --maxthreads ~{nthreads} -o somatic_~{pairName} ~{"-i " + interval}

         ## VSQR
        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_~{pairName}/~{tumPrefix}.vcf --out somatic_~{pairName}

        ## check that this is for

            if [[ -f somatic_~{pairName}/~{tumPrefix}.vcf.recal ]];
                then 
                cp somatic_~{pairName}/~{tumPrefix}.vcf.recal somatic_~{pairName}/~{tumPrefix}.recal.vcf 
                elif [[ -f somatic_~{pairName}/~{tumPrefix}.vcf ]];
                then
                cp somatic_~{pairName}/~{tumPrefix}.vcf somatic_~{pairName}/~{tumPrefix}.recal.vcf
            fi
        fi
        
        ## Run germline based calling

        if [ ~{runGerm} -eq "1" ];
        then
        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{normalBam} -CallMNVs ~{MNVcall} -gVCF false --threadbychr ~{scatterchr} --collapse true -ploidy diploid \
        --filterduplicates true --maxthreads ~{nthreads} -o somatic_~{pairName} ~{"-i " + interval}

        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_~{pairName}/~{normPrefix}.vcf --out somatic_~{pairName}

        if [[ -f somatic_~{pairName}/~{normPrefix}.vcf.recal ]];
            then
              cp somatic_~{pairName}/~{normPrefix}.vcf.recal somatic_~{pairName}/~{normPrefix}.recal.vcf
            elif [[ -f somatic_~{pairName}/~{normPrefix}.vcf ]];
            then
              cp somatic_~{pairName}/~{normPrefix}.vcf somatic_~{pairName}/~{normPrefix}.recal.vcf
        fi
        ## dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf somatic_~{pairName}/~{normPrefix}.recal.vcf --bam ~{normalBam}
        fi
        
        if [ ~{matchPair} -eq "1" ];
        then 
            dotnet /app/VennVcf_5.2.10.49/VennVcf.dll -if somatic_~{pairName}/~{tumPrefix}.recal.vcf,somatic_~{pairName}/~{normPrefix}.recal.vcf -o venn

            tar -zvcf ~{pairName}_venn_pisces.tar.gz venn

            sampID="venn/~{tumPrefix}.recal_not_~{normPrefix}.recal.vcf"

            awk '! /\#/' $sampID | awk '{if(length($4) > length($5)) print $1"\t"($2-1)"\t"($2+length($4)-1); else print $1"\t"($2-1)"\t"($2+length($5)-1)}' > output.bed

            dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{normalBam} -CallMNVs ~{MNVcall} -gVCF true --threadbychr ~{scatterchr} --collapse true -ploidy diploid --filterduplicates true --maxthreads ~{nthreads} -o variant2_~{pairName} -i output.bed

            dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf variant2_~{pairName}/~{normPrefix}.genome.vcf --out variant2_~{pairName}
            if [[ -f variant2_~{pairName}/~{normPrefix}.genome.vcf.recal ]];
            then 
            mv variant2_~{pairName}/~{normPrefix}.genome.vcf.recal variant2_~{pairName}/~{normPrefix}.genome.recal.vcf
            elif [[ -f variant2_~{pairName}/~{normPrefix}.genome.vcf ]];
            then
            cp variant2_~{pairName}/~{normPrefix}.genome.vcf variant2_~{pairName}/~{normPrefix}.genome.recal.vcf
            fi
            
            mv venn/~{tumPrefix}.recal_not_~{normPrefix}.recal.vcf ~{tumPrefix}.somatic.unique.recal.vcf
        elif [ ~{tOnly} -eq "1" ];
        then
            mv somatic_~{pairName}/~{tumPrefix}.recal.vcf ~{tumPrefix}.somatic.unique.recal.vcf
        fi; 

        if [ ~{runScylla} -eq "1" ];
        then  
        ## perform phasing here
        dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf ~{tumPrefix}.somatic.unique.recal.vcf --bam ~{tumorBam}
        fi

        if [ ~{saveDict} -eq "1" ];
        then
            tar -zvcf refPisces.tar.gz $sname
        fi
        
    >>>

    output {
        File tumor_unique_variants= select_first(["~{tumPrefix}.somatic.unique.recal.vcf", "somatic_~{pairName}/~{tumPrefix}.recal.vcf" ])
        File tumor_unique_variants_phased= "~{tumPrefix}.somatic.unique.recal.phased.vcf" 
        File? normal_variants_same_site= select_first(["variant2_~{pairName}/~{normPrefix}.genome.recal.vcf", "NULL"])
        File? venn_zip="~{pairName}_venn_pisces.tar.gz"
        File? refzip="refPisces.tar.gz"

    }

    runtime {
        docker: "trinhanne/pisces:5.2.10"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}

task runpiscesSomaticTumOnly {
    input {
    File? refFasta
    File? refFastaFai
    File? refFastaDict
    File? pisces_reference
    File? tumorBam
    File? tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    Boolean MNVcall = true
    Boolean scatterchr = true
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    String saveDict = "0"
    String runMode
    String runScylla = "0"
    }


    String normP = select_first([ normalBam ,""])
    String tumP = select_first([ tumorBam ,""])
    String tumPrefix=if (tumP!="") then basename(sub(tumP,"\\.bam$", "")) else ""
    String normPrefix= if (normP!="") then basename(sub(normP,"\\.bam$", "")) else ""
    String buildRef = if defined(pisces_reference) then "0" else "1"

    String runTum = if (runMode!="Germline") then "1" else "0"
    String runGerm = if (runMode!="TumOnly" ||  defined(normalBam)) then "1" else "0"
    String matchPair = if (runMode=="Paired" ||  ( defined(normalBam) && defined(tumorBam))) then "1" else "0"
    String tOnly = if (runMode=="TumOnly" ||  ( !defined(normalBam) && defined(tumorBam))) then "1" else "0"
    Int tumBamSize = select_first([ceil(size(tumorBam, "G")+size(tumorBai, "G")), 0])
    Int normBamSize = select_first([ceil(size(normalBam, "G")+size(normalBai, "G")), 0])
    Int disk_size=2*(tumBamSize+ normBamSize+3)
   
    command <<<
        set -e

        mkdir somatic_~{pairName}
        sname=~{pisces_reference}

        ## Build the reference libraries here if they do not exist

        if [ ~{buildRef} -eq "1" ];
        then
        echo 'create the reference'
        sname=`basename ~{refFasta}`
        sname=$(echo $sname| cut -f 1 -d '.')
        echo $sname
        mkdir $sname
        cp ~{refFasta} $sname
        cp ~{refFastaFai} $sname
        cp ~{refFastaDict} $sname
        dotnet /app/CreateGenomeSizeFile_5.2.10.49/CreateGenomeSizeFile.dll -g $sname -s "Homo sapien $sname" -o $sname
        fi 
        
        if [ ~{buildRef} -eq "0" ];
        then
        ## to unpack the tar file
        tar xvzf ~{pisces_reference}
        ## use this when using a reference file
        ## cp -r ~{pisces_reference} . 
        sname=`basename ~{pisces_reference}`
        sname=$(echo $sname| cut -f 1 -d '.')
        fi

        ## Run the tumor based calling of variants

        if [ ~{runTum} -eq "1" ];
        then        
        ## run the variant calling: this works
        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{tumorBam} -CallMNVs ~{MNVcall} -gVCF false --threadbychr ~{scatterchr} --collapse true -c 5 \
        --filterduplicates true --maxthreads ~{nthreads} -o somatic_~{pairName} ~{"-i " + interval}

         ## VSQR
        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_~{pairName}/~{tumPrefix}.vcf --out somatic_~{pairName}

        ## check that this is for

            if [[ -f somatic_~{pairName}/~{tumPrefix}.vcf.recal ]];
                then 
                cp somatic_~{pairName}/~{tumPrefix}.vcf.recal somatic_~{pairName}/~{tumPrefix}.recal.vcf 
                elif [[ -f somatic_~{pairName}/~{tumPrefix}.vcf ]];
                then
                cp somatic_~{pairName}/~{tumPrefix}.vcf somatic_~{pairName}/~{tumPrefix}.recal.vcf
            fi
        fi
        
        ## Run germline based calling

        if [ ~{runGerm} -eq "1" ];
        then
        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{normalBam} -CallMNVs ~{MNVcall} -gVCF false --threadbychr ~{scatterchr} --collapse true -ploidy diploid \
        --filterduplicates true --maxthreads ~{nthreads} -o somatic_~{pairName} ~{"-i " + interval}

        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_~{pairName}/~{normPrefix}.vcf --out somatic_~{pairName}

        if [[ -f somatic_~{pairName}/~{normPrefix}.vcf.recal ]];
            then
              cp somatic_~{pairName}/~{normPrefix}.vcf.recal somatic_~{pairName}/~{normPrefix}.recal.vcf
            elif [[ -f somatic_~{pairName}/~{normPrefix}.vcf ]];
            then
              cp somatic_~{pairName}/~{normPrefix}.vcf somatic_~{pairName}/~{normPrefix}.recal.vcf
        fi
        ## dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf somatic_~{pairName}/~{normPrefix}.recal.vcf --bam ~{normalBam}
        fi
        
        if [ ~{matchPair} -eq "1" ];
        then 
            dotnet /app/VennVcf_5.2.10.49/VennVcf.dll -if somatic_~{pairName}/~{tumPrefix}.recal.vcf,somatic_~{pairName}/~{normPrefix}.recal.vcf -o venn

            tar -zvcf ~{pairName}_venn_pisces.tar.gz venn

            sampID="venn/~{tumPrefix}.recal_not_~{normPrefix}.recal.vcf"

            awk '! /\#/' $sampID | awk '{if(length($4) > length($5)) print $1"\t"($2-1)"\t"($2+length($4)-1); else print $1"\t"($2-1)"\t"($2+length($5)-1)}' > output.bed

            dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ~{normalBam} -CallMNVs ~{MNVcall} -gVCF true --threadbychr ~{scatterchr} --collapse true -ploidy diploid --filterduplicates true --maxthreads ~{nthreads} -o variant2_~{pairName} -i output.bed

            dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf variant2_~{pairName}/~{normPrefix}.genome.vcf --out variant2_~{pairName}
            if [[ -f variant2_~{pairName}/~{normPrefix}.genome.vcf.recal ]];
            then 
            mv variant2_~{pairName}/~{normPrefix}.genome.vcf.recal variant2_~{pairName}/~{normPrefix}.genome.recal.vcf
            elif [[ -f variant2_~{pairName}/~{normPrefix}.genome.vcf ]];
            then
            cp variant2_~{pairName}/~{normPrefix}.genome.vcf variant2_~{pairName}/~{normPrefix}.genome.recal.vcf
            fi
            
            mv venn/~{tumPrefix}.recal_not_~{normPrefix}.recal.vcf ~{tumPrefix}.somatic.unique.recal.vcf
        elif [ ~{tOnly} -eq "1" ];
        then
            mv somatic_~{pairName}/~{tumPrefix}.recal.vcf ~{tumPrefix}.somatic.unique.recal.vcf
        fi; 

        if [ ~{runScylla} -eq "1" ];
        then  
        ## perform phasing here
        dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf ~{tumPrefix}.somatic.unique.recal.vcf --bam ~{tumorBam}
        fi

        if [ ~{saveDict} -eq "1" ];
        then
            tar -zvcf refPisces.tar.gz $sname
        fi
        
    >>>

    output {
        File tumor_unique_variants= select_first(["~{tumPrefix}.somatic.unique.recal.vcf", "somatic_~{pairName}/~{tumPrefix}.recal.vcf" ])
        File tumor_unique_variants_phased= "~{tumPrefix}.somatic.unique.recal.phased.vcf" 
        File? normal_variants_same_site= select_first(["variant2_~{pairName}/~{normPrefix}.genome.recal.vcf", "NULL"])
        File? venn_zip="~{pairName}_venn_pisces.tar.gz"
        File? refzip="refPisces.tar.gz"

    }

    runtime {
        docker: "trinhanne/pisces:5.2.10"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}