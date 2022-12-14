
## Script for collating information from different runs

version 1.0

import "oncokb.wdl" as oncokb

workflow ClinicalReport {
# inputs and their types specified here
    input {
        File inputSNV
        File? inputCV
        File? inputSV
        File? ploidyTar
        String sampleName
        String token
        String searchby
        String oncotree
        File cosmicMut
        File cosmicGenes
        File MsigDBAnnotation
        File pfam
        File pirsf
        File ACMG
        File? AddList
        String? grepRm
        File? CNVplot
        File inputYaml
        File GTex
        Int memoryGB = 10
        File? columnEntries
        String runMode ="Germline"
        Boolean SNVvcfformat =true
        Boolean canonical = true
        Boolean runOncokb
        Boolean CNVfunco = true
        File AAlist
        File param_config
#        Boolean runOncokb = false
        String dockerFile = "trinhanne/clin_report_annot:v2.7"
    }

        Map[String, String] filt_params=read_json(param_config)

    if (runOncokb && SNVvcfformat) {
        call oncokb.oncokb as OncokbWF {
            input:
            vcf = inputSNV,
            oncotree = oncotree,
            samplename = sampleName,
            token = token,
            searchby = searchby,
            AAlist = AAlist,
            canonical=canonical,
            runMode=runMode,
            grepRm=grepRm
        }
    }

    Boolean SNVvcfformat2 = if (defined(OncokbWF.oncokbout) || SNVvcfformat ) then false else SNVvcfformat
    File SNVin = select_first([OncokbWF.oncokbout, inputSNV])

    call ConvertSNVs {
        input:
        inputSNV = SNVin,
        sampleName = sampleName,
        cosmicMut = cosmicMut,
        cosmicGenes = cosmicGenes,
        MsigDBAnnotation = MsigDBAnnotation,
        SNVvcfformat=SNVvcfformat2,
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


    if (defined(inputCV)){
    File inputDataCV = select_first([inputCV, 'NULL'])
    call ConvertSVs as CNVFormat {
        input:
        inputSV=inputDataCV,
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
        Passfilt=filt_params["SVPass"],
        runMode=runMode,
        funco=CNVfunco
    }
    }

    if (defined(inputSV)){
    File inputDataSV = select_first([inputSV, 'NULL'])
     call ConvertSVs as SVFormat {
        input:
        inputSV=inputDataSV,
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
        Passfilt=filt_params["SVPass"],
        runMode=runMode,
        funco=false
    }}

    call CreateClinical {
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
        param_config=param_config,
        CNVplot=CNVplot,
        runMode=runMode
    }


# outputs and their types specified here
    output {
      File? yaml=CreateClinical.yamlOut
      File? CNV = CNVFormat.CNV
      File? htmlFile = select_first([CreateClinical.germline_html, CreateClinical.tumour_html])
      File FinalannotMafGz = ConvertSNVs.annotMafGz
      File SupportingSNVs = ConvertSNVs.SNV
    }

}

task CreateClinical {
    input {
        File FormatSNV
        File? CNV 
        File? SV
        File? ploidyTar
        File? SVsplit
        String sampleName
        File yaml
        Int memoryGB
        String dockerFile
        File param_config
        String runMode
        File? CNVplot
    }

    command <<<
        set -euxo pipefail
        echo 'move items and export names if they exist' 
        export runMode="~{runMode}"
        ## SNV data
        if [ -f ~{FormatSNV} ]; then
            tar -C . -xvf ~{FormatSNV}
             export snvsummary="~{sampleName}variantSummary.filt.maf"
                export snvcancer="~{sampleName}cancerVariants.filt.maf"
                export snvvus="~{sampleName}pathogenicVUS.filt.maf"
                export snvdrug="~{sampleName}drugprotective.filt.maf"
                export snvhallmark="~{sampleName}PathwayVariants.filt.maf"
                export snvacmg="~{sampleName}ACMG.filt.maf"
        else 
            export snvsummary=""
            export snvcancer=""
            export snvvus=""
            export snvdrug=""
            export snvhallmark=""
            export snvacmg=""
        fi
        ## CNV data
        if [ -f "~{CNV}" ]; then
            echo 'CNV file exists'
            cp ~{CNV} .
            export cnvAnnot="~{sampleName}.CNV.formated.tsv"
        else
            echo 'CNV file does not exists'
            export cnvAnnot=""
        fi
        ## SV data
        if [ -f "~{SV}" ]; then
            echo 'SV file exists'
            cp ~{SV} .
            export svAnnot="~{sampleName}.SV.formated.tsv"
        else
            echo 'SV file does not exist'
            export svAnnot=""
        fi
        if [ -f "~{SVsplit}" ]; then
            echo 'SV split file exists'
            cp ~{SVsplit} .
            export SVsplit="~{sampleName}.SV.split.formated.tsv"
        else 
            echo 'SV split file does not exists'
            export SVsplit=""
        fi
        ## CNV plot
        if [ -f "~{CNVplot}" ]; then
            echo 'Tumour CNV plot exists'
            cp ~{CNVplot} ~{sampleName}.modeled.png
            export CNVplot="~{sampleName}.modeled.png"
        else
            echo 'Tumour CNV plot does not exists'
            export CNVplot=""
        fi
        # copy template files over
        cp /templates/*.Rmd . 
        cp ~{param_config} ./parameter.json

        if [ ~{runMode} == "Germline" ]; then
        echo 'prepare for germline mode'
        tar -C . -xvf ~{ploidyTar}
            if [ -d SAMPLE_0 ]; then
                mv ./SAMPLE_0/contig_ploidy.tsv .
            fi
        mv Template_Germline_Report.Rmd ~{sampleName}_Germline_Report.Rmd
        export titanParams=""
        else 
        echo 'prepare for tumour mode'
        cp ~{ploidyTar} ~{sampleName}.optimalclusters.txt
        export titanParams="~{sampleName}.optimalclusters.txt"
        mv Template_Somatic_Report.Rmd ~{sampleName}_Somatic_Report.Rmd
        fi


        # export everything and edit the yaml file
        #echo $SVsplit
       ## create a yaml file for the outputs
        echo 'edit yaml' 
        rm -f final.yml temp.yml
        ( echo "cat > final.yml <<EOF"; cat "~{yaml}"; echo "EOF";)>temp.yml
        ##( echo "cat <<EOF >final.yml"; cat ~{yaml}) > temp.yml
        . temp.yml
        # print the input file
         cat final.yml
        # Now use this yaml and use it in the Rmd
        ## How to use whethere pandoc exists?>/Applications/RStudio.app/Contents/MacOS/pandoc # /usr/local/bin/pandoc
        if [ ~{runMode} == "Germline" ]; then
        Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc"); rmarkdown::render("./~{sampleName}_Germline_Report.Rmd")'   
        else 
        Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc"); rmarkdown::render("./~{sampleName}_Somatic_Report.Rmd")'   
        fi
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
        File? germline_html = "~{sampleName}_Germline_Report.html"
        File? tumour_html = "~{sampleName}_Somatic_Report.html"
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
        String runMode 
        Boolean funco
    }

    Boolean is_compressed = sub(basename(inputSV), ".*\\.", "") == "gz"
    Boolean isGermline = if (runMode=="Germline") then true else false

    command <<<
        set -euxo pipefail
        TissueWd=`grep "Tissue" ~{inputYaml}| cut -d' ' -f2 `
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `

        ## ~/gitLibs/DockerRepository/clinRep
        if [ ~{CNV} == true ] && [ ~{funco} == true ]; then
        echo 'Run tumour mode annotations funco annotated'
        Rscript /opt/AnnotateTumCNV.R --tsv ~{inputSV} --outputname ~{sampleName} --MSigDB ~{MsigDBAnnotation} \
        --GTex ~{GTex} --CosmicList ~{cosmicGenes} ~{"--AddList " + AddList} --pathwayList "$keyWd" --Tissue "$TissueWd" 
        else
        echo 'Run germline mode annotations'
        Rscript /opt/SummarizeAnnotSV.R --AnnotSVtsv ~{inputSV} --outputname ~{sampleName} --MSigDB ~{MsigDBAnnotation} \
        --GTex ~{GTex} --CosmicList ~{cosmicGenes} ~{"--AddList " + AddList} --pathwayList "$keyWd" --ACMGCutoff ~{ACMGcutoff} --Tissue "$TissueWd" --CNV ~{CNV} \
        ~{"--SRfilter " + SRcounts} ~{"--PRfilter " + PRcounts} --PASSfilt ~{Passfilt} --germline ~{isGermline}
        fi

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
        set -euxo pipefail
        annotMaf="~{sampleName}_oncokb.maf"
        tempOut="~{sampleName}_newannot.maf"
        annotM2="~{sampleName}_oncokb_final.maf"

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
