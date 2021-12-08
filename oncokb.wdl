version 1.0

# this workflow annotes a maf file with oncokb

workflow oncokbAnnotate {
# inputs and their types specified here
    input {
    File maf
    String oncotree
    String samplename
    String token
    # Options for searchby: "hgvsp", "hgvsp_short", "hgvsg"
    String searchby = "hgvsp" 
    }
# start listing all the tasks in your workflow here and the required inputs. This example only has one

    call oncokb {
        input:
        maf = maf,
        oncotree = oncotree,
        samplename = samplename,
        token = token,
        searchby = searchby 
    }
# outputs and their types specified here
    output {
        File onco = oncokb.oncokbout
    }
}



task oncokb {
    input {
    File maf
    String token
    String samplename
    String searchby
    String oncotree
    }

    command {

python <<CODE

with open('clinannot.txt', 'w') as f:
    f.write("SAMPLE_ID\tONCOTREE_CODE\n")
    f.write('${samplename}'+"\t"+'${oncotree}')
f.close()
CODE

echo `head ${maf}`

python /oncokb/MafAnnotator.py -i ${maf} -o "${samplename}_oncokb.maf" -c "clinannot.txt" -b ${token} -q ${searchby}
    }

    output {
    File oncokbout = "${samplename}_oncokb.maf"
    }

    runtime {
    docker: "trinhanne/oncokb:3.1.1"
    preemptible: "3"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    }
}