## Script for collating information from different runs

version 1.0

workflow ClinicalReportSNVsOnly {
# inputs and their types specified here
    input {
        File inputSNV
        String sampleName
        File cosmicMut
        File cosmicGenes
        File MsigDBAnnotation
        File pfam
        File pirsf
        File ACMG
        String runMode ="Tumour"
        File? AddList
        String FiltOut = "damaging|ncogen|pathogenic|risk_factor|protective|TSG|drug_response|fusion"
        File inputYaml
        Int memoryGB = 14
        File AAlist
        Boolean canonical = true
        File? columnEntries
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
        inputYaml=inputYaml,
        canonical=canonical,
        AAlist=AAlist,
        runMode=runMode,
        columnEntries=columnEntries
    }

    call CreateClinical{
        input:
        FormatSNV=ConvertSNVs.SNV,
        sampleName=sampleName,
        memoryGB=memoryGB,
        yaml=inputYaml
    }

# outputs and their types specified here
	output {
	  File yaml=CreateClinical.yamlOut
      File htmlFile = CreateClinical.html
      File FinalannotMafGz = ConvertSNVs.annotMafGz
      File SupportingSNVs = ConvertSNVs.SNV
	}

}

task CreateClinical {
    input {
        File FormatSNV
        String sampleName
        File yaml
        Int memoryGB
    }

    command <<<

        tar -C . -xvf ~{FormatSNV}
        ###tar {FormatSNV}
        cp /template/*.Rmd . 
        ##cp ~/Documents/ER_pilot/New_Clin_Reports/*.Rmd .
        mv Template_Germline_Report.Rmd ~{sampleName}_Germline_Report.Rmd

        # export everything and edit the yaml file
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
        docker: "trinhanne/clin_report_annot:v2.1"
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
        docker: "trinhanne/clin_report_annot:v2.1"
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
        String runMode
        Boolean canonical
        File inputYaml
        Int memoryGB
        Int? caddScore = 20
        Float? gnomadcutoff =0.1
        File AAlist
        File? columnEntries
    }
        Int diskGB=6*ceil(size(inputSNV, "GB")+size(cosmicMut, "GB"))


    command <<<
        ## modify the SNV files here:
        #unzip files
        annotMaf="~{sampleName}_oncokb.maf"
        tempOut="~{sampleName}_newannot.maf"
        annotM2="~{sampleName}_oncokb_final.maf"

        #gunzip -c ~{inputSNV} > $annotMaf

        ##tar -xvzf ~{inputSNV} -C .
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `

        echo 'convert vcf to maf'
        Rscript /opt/vepVCF2maf.R --vcffile ~{inputSNV} --outputfile $annotMaf --sampleName ~{sampleName} --canonical "~{canonical}" --runMode ~{runMode} --AAlist ~{AAlist} 
        echo 'MAF file saved!'
        echo 'annotating SNVs with additional databases'
        Rscript /opt/DBAnnotations.R --maffile $annotMaf --outputfile $tempOut --cosmicMut ~{cosmicMut} --MSigDB ~{MsigDBAnnotation} --pfam ~{pfam} --pirsf ~{pirsf}
        paste $annotMaf $tempOut > $annotM2
        echo 'New updated tables'
        Rscript /opt/SummarizeVariants.R --maffile $annotM2 --outputname ~{sampleName} --caddscore ~{caddScore} --gnomadcutoff ~{gnomadcutoff} \
        --pathwayList "$keyWd" --ACMG ~{ACMG} ~{"--AddList " + AddList } ~{"--ColumnEntries " + columnEntries}
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
        docker: "trinhanne/clin_report_annot:v2.1"
        preemptible: "3"
        memory: memoryGB + "GB"
        disks: "local-disk ~{diskGB} HDD"
    }

}