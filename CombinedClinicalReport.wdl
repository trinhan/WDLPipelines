## Script for collating information from different runs

version 1.0

workflow ClinicalReport {
# inputs and their types specified here
    input {
        File inputSNV
        File? inputCV
        File? inputSV
        String sampleName
        File cosmicMut
        File MsigDBAnnotation
        String FiltOut = "damaging|ncogen|pathogenic|risk_factor|protective|TSG|drug_response|fusion"
        File inputYaml
        Int memoryGB = 14
    }

# start listing all the tasks in your workflow here and the required inputs. This example only has one

    call CombineReport {
        input:
        inputSNV = inputSNV,
        inputCV = inputCV,
        inputSV = inputSV,
        sampleName = sampleName,
        cosmicMut = cosmicMut,
        MsigDBAnnotation = MsigDBAnnotation,
        FiltOut=FiltOut,
        inputYaml=inputYaml,
        memoryGB=memoryGB
    }
# outputs and their types specified here
	output {
		File yaml=CombineReport.yaml
        Array[File] SNV = CombineReport.SNV
        File htmlFile = CombineReport.htmlFile
	}

}

task CombineReport {
    input {
        File inputSNV
        File? inputCV
        File? inputSV
        String sampleName
        File cosmicMut
        File MsigDBAnnotation
        String FiltOut
        File inputYaml
        Int memoryGB
    }
        Int diskGB=6*ceil(size(inputSNV, "GB")+size(inputCV, "GB")+size(inputSV, "GB")+size(cosmicMut, "GB"))


    command <<<

        #unzip files
        annotMaf="~{sampleName}_oncokb.maf"
        tempOut="~{sampleName}_newannot.maf"
        annotM2="~{sampleName}_oncokb_2.maf"

        gunzip -c ~{inputSNV} > $annotMaf

        ##tar -xvzf ~{inputSNV} -C .
        keyWd=`grep "TreatmentKeywords" ~{inputYaml}| cut -d' ' -f2 `


        echo 'annotating SNVs with additional databases'
        Rscript /opt/DBAnnotations.R --maffile $annotMaf --outputfile $tempOut --cosmicMut ~{cosmicMut} --MSigDB ~{MsigDBAnnotation}
        paste $annotMaf $tempOut > $annotM2
        ## filter out here
        grep -E "~{FiltOut}" $annotM2> ~{sampleName}.filt.prot.onco.maf
        #echo 'summarise SNVs for export'
        Rscript /opt/SummarizeVariants.R --maffile $annotM2 --outputname ~{sampleName} --pathwayList "$keyWd"
        # compress the original file
        tar -cvzf ~{sampleName}_oncokb_2.maf.gz ~{sampleName}_oncokb_2.maf

        ## create a yaml file for the outputs

        export snvsummary="~{sampleName}variantSummary.filt.maf"
        export snvcancer="~{sampleName}cancerVariants.filt.maf"
        export snvvus="~{sampleName}pathogenicVUS.filt.maf"
        export snvdrug="~{sampleName}drugprotective.filt.maf"
        export snvhallmark="~{sampleName}HallmarkVUS.filt.maf"

        rm -f final.yml temp.yml
        ( echo "cat <<EOF >final.yml"; cat ~{inputYaml}) > temp.yml
        . temp.yml

        cp /template/*.Rmd . 
        mv Template_Germline_Report.Rmd ~{sampleName}_Germline_Report.Rmd

        # Now use this yaml and use it in the Rmd
        Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin/pandoc"); rmarkdown::render("./~{sampleName}_Germline_Report.Rmd")'

        tar -cvzf ~{sampleName}.filt.prot.onco.maf.gz ~{sampleName}.filt.prot.onco.maf

    >>>

    output {
        File yaml = "final.yml"
        Array[File] SNV = glob("*filt.maf")
        File htmlFile = "~{sampleName}_Germline_Report.html"
        File annotMafGz = "~{sampleName}_oncokb_2.maf.gz"
        File mafFiltGz = "~{sampleName}.filt.prot.onco.maf.gz"
    }

    runtime {
        docker: "trinhanne/clin_report_annot:v1"
        preemptible: "3"
        memory: memoryGB + "GB"
        disks: "local-disk ~{diskGB} HDD"
    }

}