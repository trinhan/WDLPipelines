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
    String searchby = "HGVSp_Short" 
    File pfam
    File pirsf
    File AAlist
    String? grepRm
    String? FiltOut = "damaging|ncogenic|pathogenic|risk_factor"
    String? canonical = "FALSE"
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
        grepRm=grepRm,
        FiltOut=FiltOut,
        canonical=canonical
    }
# outputs and their types specified here
    output {
        File onco = oncokb.oncokbout
        File MafFilt = oncokb.MafFilt
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
    String? FiltOut
    String? canonical = "FALSE"
    }

    Int diskGB = ceil(2*size(vcf, "G")+size(pfam ,"G"))*3 + diskGB_buffer
    String rmSamps = if defined(grepRm) then "1" else "0"


    command {

python3 <<CODE

with open('clinannot.txt', 'w') as f:
    f.write("SAMPLE_ID\tONCOTREE_CODE\n")
    f.write('${samplename}'+"\t"+'${oncotree}')
    f.close()
CODE

    OutputMaf=${samplename}.maf
    OutputMaf2=${samplename}.prot.maf
    MafFilt=${samplename}.prot.onco.filt.maf
    vcfMod=${samplename}.vcf
    vcfMod2=${samplename}.2.vcf


    gunzip -c ${vcf} > $vcfMod

    if [ ${rmSamps} -eq "1" ] ;
    then 
        grep -v -E "${grepRm}" $vcfMod > $vcfMod2
        mv $vcfMod2 $vcfMod
    fi

    ## parallelise if the number of samples is super long?

# Run the Rscript: step 1, create the maf
Rscript /opt/vepVCF2maf4Oncokb.R --vcffile $vcfMod --outputfile $OutputMaf --sampleName ${samplename} --AAlist ${AAlist} --canonical "${canonical}"
Rscript /opt/annotateProteins.R --maffile $OutputMaf --outputfile $OutputMaf2 ~{"--pfam " + pfam} ~{"--pirsf " + pirsf}
# Run the Rscript: step 2 annotate the data file with pfam and pirsf
# Run the Rscript: step 2 annotate the data file with pfam and pirsf

python3 /oncokb/MafAnnotator.py -i $OutputMaf2 -o "${samplename}_oncokb.maf" -c "clinannot.txt" -b ${token} -q ${searchby}

# Filter out 
grep -E "${FiltOut}" ${samplename}_oncokb.maf > $MafFilt

# compress 
gzip $MafFilt
gzip ${samplename}_oncokb.maf

    }

    output {
    File oncokbout = "${samplename}_oncokb.maf.gz"
    File MafFilt = "${samplename}.prot.onco.filt.maf.gz"
    }

    runtime {
    docker: "trinhanne/oncokb:3.1.2Renv"
    preemptible: "3"
    memory: memoryGB + "GB"
    disks: "local-disk ${diskGB} HDD"
    }
}