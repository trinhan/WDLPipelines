# This workflow runs pisces in a parallelised manner for Germline variant calling
# This WF takes in a list of parallelised inputs, or can generate this from scratch if inputs 

version 1.0

workflow pisces_workflow {
  input {
    File refFasta
    File refFastaIdx
    File refFastaDict
    File? pisces_reference
    File normalBam
    File normalBai
    String pairName
    File? InputtargetInterval
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
     call runpiscesGermline {
        input:
            refFasta=refFasta,
            refFastaFai=refFastaIdx,
            refFastaDict=refFastaDict,
            normalBam=normalBam,
            normalBai=normalBai,
            pairName=pairName,
            pisces_reference=pisces_reference,
            interval=bed_list[idx]
        }
    }

    call UpdateHeaders as NormHeaders {
        input:
            input_vcfs =runpiscesGermline.normal_variants,
            ref_dict=refFastaDict,
            gatk_docker = gatk_docker,
            caller = "Pisces"
    }

    call CombineVariants as NormVariants {
        input:
            input_header=NormHeaders.head_vcf,
            ref_fasta = refFasta,
            ref_fai = refFastaIdx,
            ref_dict = refFastaDict,
            gatk_docker = gatk_docker,
            sample_name = pairName,
            caller="Pisces"
    }

    output {
        File normal_variants=NormVariants.merged_vcf
    }
}

task runpiscesGermline {
    input {
    File? refFasta
    File? refFastaFai
    File? refFastaDict
    File? pisces_reference
    File normalBam
    File normalBai
    String pairName
    File? interval
    Boolean MNVcall = true
    Boolean scatterchr = true
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    String saveDict = "0"
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
        
    >>>

    output {
        File normal_variants =  select_first(["somatic_${pairName}/~{normPrefix}.recal.vcf", "somatic_~{pairName}/~{normPrefix}.vcf"])
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

