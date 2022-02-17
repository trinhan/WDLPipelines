################
## Workflow uses oncokb to annotate a vep vcf
## github: https://github.com/oncokb/oncokb-annotator. Uses v.3.1.2
## NB: This only works for SNVs right now
## ############
## Requirements:
## vcf: vep annotated vcf
## oncotree: e.g. MEL. oncotree code. See http://oncotree.mskcc.org/#/home
## token: token code required. Sign up here to get one for API access: https://www.oncokb.org/apiAccess
## searchby: options are hgvsp, hgvsp_short, hgvsg
## AAlist: File to convert 3 letter AA to 1 letter (in .json)
## grepRm (optional): To shorten the search time, remove all rows within the vcf containing these terms. Use grep terminology here. e.g. "synonymous|benign|intergenic"
## FiltOut (optional): Save a compressed maf file based on these terms .e.g. "damaging|pathogenic"
## canonical : "TRUE" or "FALSE". VEP run with multiple annotations per variant. Do you want to filter out based on canonical version?

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
    File AAlist
    String? grepRm = "intergenic|synonymous|pseudogene|non_coding_transcript_variant"
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
        grepRm=grepRm,
        canonical=canonical
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
    File AAlist
    String? memoryGB ="10"
    Float? diskGB_buffer = 5
    String? grepRm
    String? canonical = "FALSE"
    }

    Int diskGB = 3*ceil(size(vcf, "GB") + diskGB_buffer)
    Int memoryGB = 9*ceil(size(vcf, "GB"))
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

    echo ${memoryGB}
    echo ${diskGB}

    gunzip -c ${vcf} > $vcfMod

    if [ ${rmSamps} -eq "1" ] ;
    then 
        grep -v -E "${grepRm}" $vcfMod > $vcfMod2
        mv $vcfMod2 $vcfMod
    fi

    ## parallelise if the number of samples is super long?

# Run the Rscript: step 1, create the maf

    Rscript /opt/vepVCF2maf4Oncokb.R --vcffile $vcfMod --outputfile $OutputMaf --sampleName ${samplename} --canonical "${canonical}" 
    Rscript /opt/HGVSMafAnnot.R --maffile $OutputMaf --outputfile $OutputMaf2 --AAlist ${AAlist}
    ## Run the Rscript: step 2 annotate the data file with pfam and pirsf
    ## Run the Rscript: step 2 annotate the data file with pfam and pirsf

python3 /oncokb/MafAnnotator.py -i $OutputMaf2 -o "${samplename}_oncokb.maf" -c "clinannot.txt" -b ${token} -q ${searchby}
gzip ${samplename}_oncokb.maf

    }

    output {
    File oncokbout = "${samplename}_oncokb.maf.gz"

    }

    runtime {
    docker: "trinhanne/oncokb:v3.2.2Renv"
    preemptible: "3"
    memory: memoryGB + "GB"
    disks: "local-disk ${diskGB} HDD"
    }
}