version 1.0

#workflow pisces_workflow {
#  input {
#    File? refFasta
#    File? refFastaFai
#    File? pisces_reference
#    File tumorBam
#    File tumorBai
#    File? normalBam
#    File? normalBai
#    String pairName
#    File? interval
#   }
#
#     call runpisces {
#        input:
#            refFasta=refFasta,
#            refFastaFai=refFastaFai,
#            tumorBam=tumorBam,
#            normalBam=normalBam,
#            normalBai=normalBai,
#            tumorBai=tumorBai,
#            pairName=pairName,
#            pisces_reference=pisces_reference,
#            interval=interval
#        }
#
#    output {
#        File pisces_tum_phased=runpisces.tumor_variants_phased
#        File pisces_tum_recal=runpisces.tumor_variants_unique
#        File? pisces_norm_recal=runpisces.normal_variants
#        File? tumour_only=runpisces.venn_zip
#    }
#}

task runpisces {
    input {
    File? refFasta
    File? refFastaFai
    File? pisces_reference
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
    String pairName
    File? interval
    Boolean MNVcall = true
    Boolean scatterchr = true
    Int? nthreads =2
    String mem =8
    Int preemptible =3
    }

    Int disk_size=2*ceil(size(tumorBam, "GB")+3)
    String normP = select_first([ normalBam ,""])

    String tumPrefix=basename(sub(tumorBam,"\\.bam$", ""))
    String normPrefix= if (normP!="") then basename(sub(normP,"\\.bam$", "")) else ""
   
    command {
        set -e

        mkdir somatic_${pairName}
        echo ${normPrefix}
        sname=${pisces_reference}

        if [[ -f "${refFasta}" ]];
        then
        echo 'create the reference'
        sname=`basename ${refFasta}`
        sname=$(echo $sname| cut -f 1 -d '.')
        echo $sname
        mkdir $sname
        mv ${refFasta} $sname
        mv ${refFastaFai} $sname
        dotnet /app/CreateGenomeSizeFile_5.2.10.49/CreateGenomeSizeFile.dll -g $sname -s "Homo sapien $sname" -o $sname
        fi 

        ## run the variant calling: this works
        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ${tumorBam} -CallMNVs ${MNVcall} -gVCF false --threadbychr ${scatterchr} --collapse true -c 5 \
        --filterduplicates true --maxthreads ${nthreads} -o somatic_${pairName} ~{"-i " + interval}

         ## VSQR
        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_${pairName}/${tumPrefix}.vcf --out somatic_${pairName}

        if [[ -f somatic_${pairName}/${tumPrefix}.vcf.recal ]];
            then 
            mv somatic_${pairName}/${tumPrefix}.vcf.recal somatic_${pairName}/${tumPrefix}.recal.vcf 
            else 
            cp somatic_${pairName}/${tumPrefix}.vcf somatic_${pairName}/${tumPrefix}.recal.vcf
        fi
        
        ## phasing
        ## dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf somatic_${pairName}/${tumPrefix}.recal.vcf --bam ${tumorBam} --filterduplicates true
        ## If germline is present, do an intersection 

        if [[ -f ${normalBam} ]];
        then
        dotnet /app/Pisces_5.2.10.49/Pisces.dll -g $sname -b ${normalBam} -CallMNVs ${MNVcall} -gVCF false --threadbychr ${scatterchr} --collapse true -ploidy diploid \
        --filterduplicates true --maxthreads ${nthreads} -o somatic_${pairName} ~{"-i " + interval}

        dotnet /app/VariantQualityRecalibration_5.2.10.49/VariantQualityRecalibration.dll --vcf somatic_${pairName}/${normPrefix}.vcf --out somatic_${pairName}

        if [[ -f somatic_${pairName}/${normPrefix}.vcf.recal ]];
            then
              mv somatic_${pairName}/${normPrefix}.vcf.recal somatic_${pairName}/${normPrefix}.recal.vcf
        else
            cp somatic_${pairName}/${normPrefix}.vcf somatic_${pairName}/${normPrefix}.recal.vcf
        fi
        
        ## dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf somatic_${pairName}/${normPrefix}.recal.vcf --bam ${normalBam}
        fi
        
        if [[ -f ${normalBam} ]];
        then 
            dotnet /app/VennVcf_5.2.10.49/VennVcf.dll -if somatic_${pairName}/${tumPrefix}.recal.phased.vcf,somatic_${pairName}/${normPrefix}.recal.vcf -o venn
            tar -zvcf ~{pairName}_venn_pisces.tar.gz venn
            mv venn/${tumPrefix}.recal.phased_not_${normPrefix}.recal.vcf ${tumPrefix}.somatic.unique.recal.vcf
        else 
            mv somatic_${pairName}/${tumPrefix}.recal.vcf ${tumPrefix}.somatic.unique.recal.vcf
        fi 

        # perform phasing here
        dotnet /app/Scylla_5.2.10.49/Scylla.dll -g $sname --vcf ${tumPrefix}.somatic.unique.recal.vcf --bam ${tumorBam}
    }

    output {
        File tumor_variants_unique="${tumPrefix}.somatic.unique.recal.vcf"
        File tumor_variants_phased="${tumPrefix}.somatic.unique.recal.phased.vcf"
        File? normal_variants = "somatic_${pairName}/${normPrefix}.recal.vcf"
        File? venn_zip="~{pairName}_venn_pisces.tar.gz"
    }

    runtime {
        docker: "trinhanne/pisces:5.2.10"
        cpu: nthreads
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
    }

}

