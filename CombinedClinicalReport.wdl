
## Script for collating information from different runs

version 1.0

workflow ClinicalReport {
# inputs and their types specified here
    input {
        File inputSNV
        File inputCV
        File inputSV
        File ploidyTar
        String sampleName
        File cosmicMut
        File cosmicGenes
        File MsigDBAnnotation
        File pfam
        File pirsf
        File ACMG
        File? AddList
        File inputYaml
        File GTex
        Int memoryGB = 10
        File? columnEntries
        String runMode ="Germline"
        Boolean SNVvcfformat =true
        Boolean canonical = true
        File? AAlist
        File param_config
        String dockerFile = "trinhanne/clin_report_annot:v2.5"
    }

        Map[String, String] filt_params=read_json(param_config)

    call ConvertSNVs {
        input:
        inputSNV = inputSNV,
        sampleName = sampleName,
        cosmicMut = cosmicMut,
        cosmicGenes = cosmicGenes,
        MsigDBAnnotation = MsigDBAnnotation,
        SNVvcfformat=SNVvcfformat,
        pfam=pfam,
        pirsf=pirsf,
        memoryGB=memoryGB,
        ACMG=ACMG,
        AddList=AddList,
        inputYaml=inputYaml,
        columnEntries=columnEntries,
        canonical=canonical,
        AAlist=AAlist,
        runMode=runMode,
        dockerFile=dockerFile,
        caddScore=filt_params["OverallCADD"],
        gnomadcutoff=filt_params["OverallGnomad"],
        Agnomad=filt_params["TierAGnomad"],
        Bgnomad=filt_params["TierBGnomad"],
        Dgnomad=filt_params["TierDGnomad"],
        Bcadd=filt_params["TierBCadd"],
        Dcadd=filt_params["TierDCadd"],
        AOnlyCoding=filt_params["TierAOnlyCoding"],
        BOnlyCoding=filt_params["TierBOnlyCoding"],
        DOnlyCoding=filt_params["TierDOnlyCoding"],
        BPathogenic=filt_params["TierBPathogenic"],
        DPathogenic=filt_params["TierDPathogenic"]

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
        ACMGcutoff=filt_params["AnnotSVACMGcutoff"],
        CNV=true,
        dockerFile=dockerFile,
        Passfilt=filt_params["SVPass"]
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
        ACMGcutoff=filt_params["AnnotSVACMGcutoff"],
        CNV=false,
        dockerFile=dockerFile,
        SRcounts=filt_params["SVSRfilter"],
        PRcounts=filt_params["SVPRfilter"],
        Passfilt=filt_params["SVPass"]
    }

    call CreateClinical{
        input:
        FormatSNV=ConvertSNVs.SNV,
        sampleName=sampleName,
        CNV=CNVFormat.CNV,
        SV=SVFormat.SV,
        SVsplit=SVFormat.SVsplit,
        memoryGB=memoryGB,
        yaml=inputYaml,
        dockerFile=dockerFile,
        ploidyTar=ploidyTar,
        param_config=param_config
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
        File ploidyTar
        File SVsplit
        String sampleName
        File yaml
        Int memoryGB
        String dockerFile
        File param_config
    }

    command <<<
        echo 'move items' 
        tar -C . -xvf ~{FormatSNV}
        tar -C . -xvf ~{ploidyTar}
        ###tar {FormatSNV}
        cp ~{CNV} .
        cp ~{SV} .
        cp ~{SVsplit} .
        cp /template/*.Rmd . 
        cp ~{param_config} ./parameter.json
        ##cp ~/Documents/ER_pilot/New_Clin_Reports/*.Rmd .
        mv Template_Germline_Report.Rmd ~{sampleName}_Germline_Report.Rmd
        echo 'export items and edit yaml' 
        # export everything and edit the yaml file
        export cnvAnnot="~{sampleName}.CNV.formated.tsv"
        export svAnnot="~{sampleName}.SV.formated.tsv"
        export SVsplit="~{sampleName}.SV.split.formated.tsv"
        #echo $SVsplit
       ## create a yaml file for the outputs
        export snvsummary="~{sampleName}variantSummary.filt.maf"
        export snvcancer="~{sampleName}cancerVariants.filt.maf"
        export snvvus="~{sampleName}pathogenicVUS.filt.maf"
        export snvdrug="~{sampleName}drugprotective.filt.maf"
        export snvhallmark="~{sampleName}PathwayVariants.filt.maf"
        export snvacmg="~{sampleName}ACMG.filt.maf"
        echo 'edit yaml' 
        rm -f final.yml temp.yml
        ( echo "cat <<EOF >final.yml"; cat ~{yaml}) > temp.yml
        . temp.yml
        # print the input file
        cat final.yml
        # Now use this yaml and use it in the Rmd
        ## How to use whethere pandoc exists?>/Applications/RStudio.app/Contents/MacOS/pandoc # /usr/local/bin/pandoc
        Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc"); rmarkdown::render("./~{sampleName}_Germline_Report.Rmd")'   
        mv final.yml ~{sampleName}.final.yml
    >>>

    runtime {
        docker: dockerFile
        preemptible: "2"
        memory: memoryGB + "GB"
        disks: "local-disk 10 HDD"
        maxRetries: "1"
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
        Int? SRcounts
        Int? PRcounts 
        String? Passfilt
        String dockerFile
    }

    Boolean is_compressed = sub(basename(inputSV), ".*\\.", "") == "gz"

    command <<<

        TissueWd=`grep "Tissue" ~{inputYaml}| cut -d' ' -f2 `
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `

        ## ~/gitLibs/DockerRepository/clinRep
        Rscript /opt/SummarizeAnnotSV.R --AnnotSVtsv ~{inputSV} --outputname ~{sampleName} --MSigDB ~{MsigDBAnnotation} \
        --GTex ~{GTex} --CosmicList ~{cosmicGenes} ~{"--AddList " + AddList} --pathwayList "$keyWd" --ACMGCutoff ~{ACMGcutoff} --Tissue "$TissueWd" --CNV ~{CNV} \
        ~{"--SRfilter " + SRcounts} ~{"--PRfilter " + PRcounts} --PASSfilt ~{Passfilt}

        if [ ~{CNV} == true ]; then
            touch "~{sampleName}.SV.formated.tsv"
            touch "~{sampleName}.SV.split.formated.tsv"
        else
            touch "~{sampleName}.CNV.formated.tsv"
        fi
    >>>

    runtime {
        docker: dockerFile
        preemptible: "2"
        memory: memoryGB + "GB"
        disks: "local-disk 10 HDD"
    }

    output {
        File SV = "~{sampleName}.SV.formated.tsv"
        File SVsplit = "~{sampleName}.SV.split.formated.tsv"
        File CNV = "~{sampleName}.CNV.formated.tsv"
    }
}

task ConvertSNVs {
    input {
        File inputSNV
        String sampleName
        File cosmicMut
        File cosmicGenes
        File MsigDBAnnotation
        File ACMG
        File pfam
        File pirsf
        File? AddList
        File inputYaml
        Int memoryGB
        Int? caddScore
        Float? gnomadcutoff
        File? columnEntries
        File? pathwayFile
        String runMode
        Boolean SNVvcfformat
        File? AAlist
        Boolean canonical
        String dockerFile
        Float? Agnomad
        Float? Bgnomad
        Float? Dgnomad
        Int? Bcadd
        Int? Dcadd
        String? AOnlyCoding
        String? BOnlyCoding
        String? DOnlyCoding
        String? BPathogenic
        String? DPathogenic

    }
        Int diskGB=8*ceil(size(inputSNV, "GB")+size(cosmicMut, "GB")+size(cosmicGenes, "GB")+size(MsigDBAnnotation, "GB"))
        String AAb = select_first([AAlist, "/annotFiles/AminoAcid_table.csv"])
        String CNs = select_first([columnEntries, "/annotFiles/ColumnIDs.csv"])
        String PWF = select_first([pathwayFile, "/annotFiles/PathwayList.csv"])

    command <<<
        ## modify the SNV files here:
        #unzip files

        annotMaf="~{sampleName}_oncokb.maf"
        tempOut="~{sampleName}_newannot.maf"
        annotM2="~{sampleName}_oncokb_final.maf"

        ext="${inputSNV#*.}"
        echo $ext 

        if [ ~{SNVvcfformat} == true ]; then
            echo 'convert vcf to maf'
            Rscript /opt/vepVCF2maf.R --vcffile ~{inputSNV} --outputfile $annotMaf --sampleName ~{sampleName} --canonical ~{canonical} --runMode ~{runMode} --AAlist ~{AAb} 
            echo 'MAF file saved!'
        else     
            gunzip -c ~{inputSNV} > $annotMaf
        fi 

        ##tar -xvzf ~{inputSNV} -C .
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `

        echo 'annotating SNVs with additional databases'
        Rscript /opt/DBAnnotations.R --maffile $annotMaf --outputfile $tempOut --cosmicMut ~{cosmicMut} --cosmicGenes ~{cosmicGenes} --MSigDB ~{MsigDBAnnotation} --pfam ~{pfam} --pirsf ~{pirsf}
        paste $annotMaf $tempOut > $annotM2
        echo 'create filtered lists'
        Rscript /opt/SummarizeVariants.R --maffile $annotM2 --outputname ~{sampleName} ~{"--caddscore " + caddScore} ~{"--gnomadcutoff " + gnomadcutoff} \
        ~{"--AddList " + AddList } ~{"--columnEntries " + CNs}
        echo 'run ACMG filter'
        Rscript /opt/FindACMG.R --maffile ~{sampleName}variantsAll.maf  --outputname ~{sampleName} --ACMG ~{ACMG}
        echo 'run COSMIC filter'
        Rscript /opt/FindCosmic.R --maffile ~{sampleName}variantsAll.maf --outputname ~{sampleName} ~{"--gnomadcutoff " + Agnomad} --onlyCoding ~{AOnlyCoding}
        echo 'run Pathway filter'
        Rscript /opt/FindPathway.R --maffile ~{sampleName}variantsAll.maf  --outputname ~{sampleName} --pathwayFile ~{PWF} --pathwayList "$keyWd" \
        ~{"--gnomadcutoff " + Bgnomad} --onlyCoding ~{BOnlyCoding} ~{"--caddscore " + Bcadd} --pathogenic ~{BPathogenic}
        echo 'run Pathway filter'
        Rscript /opt/FindDrugResponse.R --maffile ~{sampleName}variantsAll.maf --outputname ~{sampleName}
        echo 'run VUS filter'
        Rscript /opt/FindVUS.R --maffile ~{sampleName}variantsAll.maf --outputname ~{sampleName} --ACMG ~{ACMG} --pathwayFile ~{PWF} --pathwayList "$keyWd" \
         ~{"--gnomadcutoff " + Dgnomad} --onlyCoding ~{DOnlyCoding} ~{"--caddscore " + Dcadd} --pathogenic ~{DPathogenic}

        tar -cvzf ~{sampleName}_oncokb_final.maf.gz ~{sampleName}_oncokb_final.maf
        ###mv *filt.maf filtOut
        tar -cvf ~{sampleName}.SNV.tar.gz *filt.maf
        ##tar -zcvf ~{sampleName}.SNV.tar.gz filtOut

    >>>

    output {
        File SNV = "~{sampleName}.SNV.tar.gz"
        File annotMafGz = "~{sampleName}_oncokb_final.maf.gz"
    }

    runtime {
        docker: dockerFile
        preemptible: "2"
        memory: memoryGB + "GB"
        disks: "local-disk ~{diskGB} HDD"
    }

}