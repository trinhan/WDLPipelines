## Script for collating information from different runs

version 1.0

workflow ClinicalReport {
# inputs and their types specified here
    input {
        File inputSNV
        File inputCV
        File inputSV
        String sampleName
        File cosmicMut
        File cosmicGenes
        File MsigDBAnnotation
        File pfam
        File pirsf
        File ACMG
        Int SVACMGcutoff =3
        File? AddList
        String FiltOut = "damaging|ncogen|pathogenic|risk_factor|protective|TSG|drug_response|fusion"
        File inputYaml
        File GTex
        Int memoryGB = 14
    }

    call ConvertSNVs {
        input:
        inputSNV = inputSNV,
        sampleName = sampleName,
        cosmicMut = cosmicMut,
        MsigDBAnnotation = MsigDBAnnotation,
        pfam=pfam,
        pirsf=pirsf,
        FiltOut=FiltOut,
        memoryGB=memoryGB,
        ACMG=ACMG,
        AddList=AddList,
        inputYaml=inputYaml
    }


    call ConvertSVs as CNVFormat {
        input:
        inputSV=inputCV,
        sampleName=sampleName,
        AddList=AddList,
        MsigDBAnnotation=MsigDBAnnotation,
        GTex=GTex,
        inputYaml=inputYaml,
        memoryGB=memoryGB,
        cosmicGenes=cosmicGenes,
        ACMGcutoff=SVACMGcutoff,
        CNV=true
    }

     call ConvertSVs as SVFormat {
        input:
        inputSV=inputSV,
        sampleName=sampleName,
        AddList=AddList,
        MsigDBAnnotation=MsigDBAnnotation,
        GTex=GTex,
        inputYaml=inputYaml,
        memoryGB=memoryGB,
        cosmicGenes=cosmicGenes,
        ACMGcutoff=SVACMGcutoff,
        CNV=false
    }

    call CreateClinical{
        input:
        FormatSNV=ConvertSNVs.SNV,
        sampleName=sampleName,
        CNV=CNVFormat.CNV,
        SV=SVFormat.SV,
        memoryGB=memoryGB,
        yaml=inputYaml
    }

# outputs and their types specified here
	output {
	  File yaml=CreateClinical.yamlOut
      File CNV = CNVFormat.CNV
      File htmlFile = CreateClinical.html
      File FinalannotMafGz = ConvertSNVs.annotMafGz
      File SupportingSNVs = ConvertSNVs.SNV
	}

}

task CreateClinical {
    input {
        File FormatSNV
        File CNV 
        File SV
        String sampleName
        File yaml
        Int memoryGB
    }

    command <<<

        tar -C . -xvf ~{FormatSNV}
        ###tar {FormatSNV}
        cp ~{CNV} .
        cp ~{SV} .
        cp /template/*.Rmd . 
        ##cp ~/Documents/ER_pilot/New_Clin_Reports/*.Rmd .
        mv Template_Germline_Report.Rmd ~{sampleName}_Germline_Report.Rmd

        # export everything and edit the yaml file
        export cnvAnnot="~{sampleName}.CNV.formated.tsv"
        export svAnnot="~{sampleName}.SV.formated.tsv"
       ## create a yaml file for the outputs
        export snvsummary="~{sampleName}variantSummary.filt.maf"
        export snvcancer="~{sampleName}cancerVariants.filt.maf"
        export snvvus="~{sampleName}pathogenicVUS.filt.maf"
        export snvdrug="~{sampleName}drugprotective.filt.maf"
        export snvhallmark="~{sampleName}PathwayVariants.filt.maf"
        export snvacmg="~{sampleName}ACMG.filt.maf"



        rm -f final.yml temp.yml
        ( echo "cat <<EOF >final.yml"; cat ~{yaml}) > temp.yml
        . temp.yml


        # Now use this yaml and use it in the Rmd
        ## How to use whethere pandoc exists?>/Applications/RStudio.app/Contents/MacOS/pandoc # /usr/local/bin/pandoc
        Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin/pandoc"); rmarkdown::render("./~{sampleName}_Germline_Report.Rmd")'
   
        mv final.yml ~{sampleName}.final.yml
    >>>

    runtime {
        docker: "trinhanne/clin_report_annot:v1.1"
        preemptible: "3"
        memory: memoryGB + "GB"
        disks: "local-disk 10 HDD"
        maxRetries: "3"
    }

    output {
        File html = "~{sampleName}_Germline_Report.html"
        File yamlOut = "~{sampleName}.final.yml"
    }

}

task ConvertSVs {
    input {
        File inputSV
        String sampleName
        File cosmicGenes
        File MsigDBAnnotation
        Int ACMGcutoff
        File GTex
        File inputYaml
        File? AddList
        Boolean CNV
        Int memoryGB 
    }

    command <<<

        TissueWd=`grep "Tissue" ~{inputYaml}| cut -d' ' -f2 `
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `

        ## ~/gitLibs/DockerRepository/clinRep
        Rscript /opt/SummarizeAnnotSV.R --AnnotSVtsv ~{inputSV} --outputname ~{sampleName} --MSigDB ~{MsigDBAnnotation} \
        --GTex ~{GTex} --CosmicList ~{cosmicGenes} ~{"--AddList " + AddList} --pathwayList "$keyWd" --ACMGCutoff ~{ACMGcutoff} --Tissue "$TissueWd" --CNV ~{CNV} --PASSfilt TRUE

        if [ ~{CNV} == true ]; then
            touch "~{sampleName}.SV.formated.tsv"
        else
            touch "~{sampleName}.CNV.formated.tsv"
        fi
    >>>

    runtime {
        docker: "trinhanne/clin_report_annot:v1.1"
        preemptible: "3"
        memory: memoryGB + "GB"
        disks: "local-disk 10 HDD"
    }

    output {
        File SV = "~{sampleName}.SV.formated.tsv"
        File CNV = "~{sampleName}.CNV.formated.tsv"
    }
}

task ConvertSNVs {
    input {
        File inputSNV
        String sampleName
        File cosmicMut
        File MsigDBAnnotation
        File ACMG
        File pfam
        File pirsf
        File? AddList
        String FiltOut
        File inputYaml
        Int memoryGB
        Int? caddScore = 20
        Float? gnomadcutoff =0.1
    }
        Int diskGB=6*ceil(size(inputSNV, "GB")+size(cosmicMut, "GB"))


    command <<<
        ## modify the SNV files here:
        #unzip files
        annotMaf="~{sampleName}_oncokb.maf"
        tempOut="~{sampleName}_newannot.maf"
        annotM2="~{sampleName}_oncokb_final.maf"

        gunzip -c ~{inputSNV} > $annotMaf

        ##tar -xvzf ~{inputSNV} -C .
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `

        echo 'annotating SNVs with additional databases'
        Rscript /opt/DBAnnotations.R --maffile $annotMaf --outputfile $tempOut --cosmicMut ~{cosmicMut} --MSigDB ~{MsigDBAnnotation} --pfam ~{pfam} --pirsf ~{pirsf}
        paste $annotMaf $tempOut > $annotM2

        Rscript /opt/SummarizeVariants.R --maffile $annotM2 --outputname ~{sampleName} --caddscore ~{caddScore} --gnomadcutoff ~{gnomadcutoff} \
        --pathwayList "$keyWd" --ACMG ~{ACMG} ~{"--AddList " + AddList }
        # compress the original file
        tar -cvzf ~{sampleName}_oncokb_final.maf.gz ~{sampleName}_oncokb_final.maf
        ###mv *filt.maf filtOut
        tar -cvf ~{sampleName}.SNV.tar.gz *.filt.maf
        ##tar -zcvf ~{sampleName}.SNV.tar.gz filtOut

    >>>

    output {
        File SNV = "~{sampleName}.SNV.tar.gz"
        File annotMafGz = "~{sampleName}_oncokb_final.maf.gz"
    }

    runtime {
        docker: "trinhanne/clin_report_annot:v1.1"
        preemptible: "3"
        memory: memoryGB + "GB"
        disks: "local-disk ~{diskGB} HDD"
    }

}