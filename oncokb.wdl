version 1.0

# this workflow annotes a maf file with oncokb

workflow oncokbAnnotate {
# inputs and their types specified here
    input {
    File vcf
    String oncotree
    String samplename
    String token
    # Options for searchby: "hgvsp", "hgvsp_short", "hgvsg"
    String searchby = "hgvsp_short" 
    File pfam
    File pirsf
    File AAlist
    String? grepRm
    }
# start listing all the tasks in your workflow here and the required inputs. This example only has one

    call oncokb {
        input:
        vcf = vcf,
        oncotree = oncotree,
        samplename = samplename,
        token = token,
        searchby = searchby,
        AAlist = AAlist,
        pfam=pfam,
        pirsf=pirsf,
        grepRm=grepRm
    }
# outputs and their types specified here
    output {
        File onco = oncokb.oncokbout
    }
}

task oncokb {
    input {
    File vcf
    String token
    String samplename
    String searchby
    String oncotree
    File pfam
    File pirsf
    File AAlist
    String? memoryGB ="10"
    String? diskGB_buffer = "20"
    String? grepRm
    }

    Int diskGB = ceil(2*size(vcf, "G")+size(pfam ,"G")+size(pirsf, "G"))*3 + diskGB_buffer
    String rmSamps = if defined(grepRm) then "1" else "0"


    command {

python3 <<CODE

with open('clinannot.txt', 'w') as f:
    f.write("SAMPLE_ID\tONCOTREE_CODE\n")
    f.write('${samplename}'+"\t"+'${oncotree}')
    f.close()
CODE

    OutputMaf=${samplename}.maf

    vcfMod=${samplename}.vcf
    vcfMod2=${samplename}.2.vcf

    gunzip -c ${vcf} > $vcfMod

    if [ ${rmSamps} -eq "1" ] ;
    then 
        grep -v -E "${grepRm}" $vcfMod > $vcfMod2
        mv $vcfMod2 $vcfMod
    fi

    ## parallelise if the number of samples is super long?
    nx2 = "nrows.out"
    `cat $vcfMod | wc -l` > $nx2
    echo $nx2

# Run the Rscript

Rscript /opt/vepVCF2maf4Oncokb.R --vcffile $vcfMod --outputfile $OutputMaf --sampleName ${samplename} --AAlist ${AAlist} --protein TRUE --pfam ${pfam} --pirsf ${pirsf}

python3 /oncokb/MafAnnotator.py -i $OutputMaf -o "${samplename}_oncokb.maf" -c "clinannot.txt" -b ${token} -q ${searchby}

    }

    output {
    File oncokbout = "${samplename}_oncokb.maf"
    }

    runtime {
    docker: "trinhanne/oncokb:3.1.2Renv"
    preemptible: "3"
    memory: memoryGB + "GB"
    disks: "local-disk ${diskGB} HDD"
    }
}