### AnnotSV: 3.1: takes a SV output and annotates this. Does not use information from GeneCards
### Requirements:
### genome Build - GRCh37 or GRCh38
### annotation file: saved in kccg-cb-tools
### snvIndelPass: To only use variants from VCF input files that passed all filters during the calling (FILTER column value equal to PASS). Default 1
### sampleName
### caller - if the caller wants to be included inthe file name
### typeofAnnotation - "full" [default] "both" or "split" [each gene takes 1 line]

version 1.0


workflow RunAnnotSV{
  input {
    File input_vcf
    File? input_vcf_idx
    String genome_build ## GRCh37,GRCh38
    File annotSVtar 
    Array[File]? snps_vcf
    String snvIndelPASS = "1"
    String sampleName
    String? caller
    String typeofAnnotation ="full"
    Int space_buffer=10
    Int memory
    String? additional_params 
  }

  call annotsv {
    input:
     input_vcf=input_vcf,
     input_vcf_idx=input_vcf_idx,
     genome_build=genome_build,
     annotSVtar=annotSVtar,
     snps_vcf=snps_vcf,
     snvIndelPASS=snvIndelPASS,
     sampleName=sampleName,
     typeofAnnotation=typeofAnnotation,
     caller=caller,
     space_buffer=space_buffer,
     memory=memory,
     additional_params=additional_params
  }

  output {
    File tsv_anotation = annotsv.sv_variants_tsv
  }
}


task annotsv {
  input {
    String genome_build
    File input_vcf
    File? input_vcf_idx
    String sampleName 
    String? caller
    Array[File]? snps_vcf
    String snvIndelPASS ="1"
    File annotSVtar
    String sampleName
    String typeofAnnotation ="full"
    Int space_buffer = 10
    Int memory = 8
    String? additional_params
  }

  Int space_needed_gb = space_buffer + ceil( size(input_vcf, "GB")+ size(annotSVtar, "GB"))

  runtime {
    memory: "~{memory} GB"
    docker: "trinhanne/annotsv:3.1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<

    # untar the file
   mkdir AnnotationsFolder
   tar xvzf ~{annotSVtar} -C AnnotationsFolder

    /opt/AnnotSV_3.1/bin/AnnotSV -SVinputFile ~{input_vcf} -bedtools /opt/bedtools2/bin/bedtools -bcftools /opt/bcftools-1.13/bcftools -snvIndelPASS ~{snvIndelPASS} \
    -genomeBuild ~{genome_build} -annotationsDir AnnotationsFolder -annotationMode ~{typeofAnnotation} -outputFile ~{sampleName}.~{caller}.annotSV.tsv -outputDir . ~{additional_params}

    gzip ~{sampleName}.~{caller}.annotSV.tsv
    #### -vcfFiles ~{sep="," snps_vcf}
  >>>

  output {
    File sv_variants_tsv = "~{sampleName}.~{caller}.annotSV.tsv.gz"
  }
}