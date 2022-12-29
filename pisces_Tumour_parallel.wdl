version 1.0

workflow pisces_workflow {
  input {
    File refFasta
    File refFastaIdx
    File refFastaDict
    File? pisces_reference
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    File? InputtargetInterval
    String runMode
    Array[File]? bed_list_in
    Array[Int]? scatterIndices_in 
    String gatk_docker
   }


   Boolean buildIndices = if defined(bed_list_in) then false else true
   File targetIntervals = select_first([InputtargetInterval, "NULL"])

if (buildIndices){
    call CallSomaticMutations_Prepare_Task {
        input:
            refFasta=refFasta,
            refFastaIdx=refFastaIdx,
            refFastaDict=refFastaDict,
            targetIntervals=targetIntervals,
            gatk_docker=gatk_docker # takes padded interval file (10bp on each side)
    }
}

    Array[File] bed_list=select_first([bed_list_in, CallSomaticMutations_Prepare_Task.bed_list])
    Array[Int] scatterIndices=select_first([scatterIndices_in, CallSomaticMutations_Prepare_Task.scatterIndices])

    scatter (idx in scatterIndices){
     call runpiscesSomatic {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaIdx,
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

    call UpdateHeaders {
        input:
            input_vcfs = runpiscesSomatic.tumor_unique_variants,
            ref_dict=refFastaDict,
            gatk_docker = gatk_docker,
            caller = "Pisces_Tum"
    }

    call CombineVariants {
        input:
            input_header=UpdateHeaders.head_vcf,
            ref_fasta = refFasta,
            ref_fai = refFastaIdx,
            ref_dict = refFastaDict,
            gatk_docker = gatk_docker,
            sample_name = pairName,
            caller="Pisces"
    }


    output {
        File tumor_variants=CombineVariants.merged_vcf
    }
}

task runpiscesSomatic {
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

        
    >>>

    output {
        File tumor_unique_variants= select_first(["~{tumPrefix}.somatic.unique.recal.vcf", "somatic_~{pairName}/~{tumPrefix}.recal.vcf" ])
    }

    runtime {
        docker: "trinhanne/pisces:5.2.10"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}

task CombineVariants {
    input {
        Array[File] input_header
        String caller
        File ref_fasta
        File ref_fai
        File ref_dict
        # runtime
        String gatk_docker
        String sample_name
        Int mem_gb= 6
    }
        Int diskGB = 4*ceil(size(ref_fasta, "GB")+size(input_header, "GB"))

    command <<<
       gatk GatherVcfs -I ~{sep=' -I ' input_header} -R ~{ref_fasta} -O ~{sample_name}.~{caller}.vcf
    >>>

    runtime {
        docker: gatk_docker
        memory: "~{mem_gb} GB"
        disk_space: "local-disk ~{diskGB} HDD"
    }

    output {
    File merged_vcf = "~{sample_name}.~{caller}.vcf"
    }
}

task UpdateHeaders {
    input {
        Array[File] input_vcfs
        File ref_dict
        # runtime
        String gatk_docker
        String caller
        Int mem_gb=6
    }
        Int diskGB = 4*ceil(size(ref_dict, "GB"))

    command <<<

        count=0
        for i in ~{sep=' ' input_vcfs}; 
        do 
         newstr=`basename $i`
         gatk UpdateVCFSequenceDictionary \
            -V $i \
            --source-dictionary ~{ref_dict} \
            --output $newstr.$count.reheader.~{caller}.vcf \
            --replace true
         count+=1
        done

    >>>

    runtime {
        docker: gatk_docker
        memory: "~{mem_gb} GB"
        disk_space: "local-disk ~{diskGB} HDD"
    }

    output {
    Array[File] head_vcf = glob("*.reheader.~{caller}.vcf")
    }
}

task CallSomaticMutations_Prepare_Task {
    input {
    # TASK INPUT PARAMS
    File targetIntervals
    File refFasta
    File refFastaIdx
    File refFastaDict

    String nWay = "10"

    # RUNTIME INPUT PARAMS
    String preemptible = "1"
    String diskGB_boot = "15"

    String gatk_docker
}

    parameter_meta {
        nWay : "Number of ways to scatter (MuTect1 and MuTect2)"
        targetIntervals : "a list of genomic intervals over which MuTect1 will operate"
        refFasta : "FASTA file for the appropriate genome build (Reference sequence file)"
        refFastaIdx : "FASTA file index for the reference genome"
        refFastaDict : "FASTA file dictionary for the reference genome"
    }

    command {
        set -euxo pipefail
        seq 0 $((~{nWay}-1)) > indices.dat
        # create a list of intervalfiles
        mkdir intervalfolder
        gatk SplitIntervals -R ~{refFasta} -L ~{targetIntervals} --scatter-count ~{nWay} --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION -O intervalfolder
        cp intervalfolder/*.interval_list .

        ## make the list of bed files
        mkdir bedfolder
        for file in *.interval_list;
        do 
            gatk IntervalListToBed -I $file -O bedfolder/$file.bed
            ##small hack to subtract 1 from the bed file
        done 

        cp bedfolder/*.bed .
    }

    runtime {
        docker         : gatk_docker
        bootDiskSizeGb : diskGB_boot
        preemptible    : preemptible
        memory         : "1 GB"
    }

    output {
        Array[File] interval_files=glob("*.interval_list")
        Array[Int] scatterIndices=read_lines("indices.dat")
        Array[File] bed_list=glob("*.bed")
    }
}




