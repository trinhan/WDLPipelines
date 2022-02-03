### AnnotSV: takes a SV output and annotates this. Does not use information from GeneCards
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
    File input_vcf_idx
    String genome_build ## GRCh37,GRCh38
    File annotSVtar 
    Array[File]? snps_vcf
    String snvIndelPASS = "1"
    String sampleName
    String? caller
    String typeofAnnotation ="full"
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
     typeofAnnotation=typeofAnnotation
  }

  output {
    File tsv_anotation = annotsv.sv_variants_tsv
  }

}


task annotsv {
  input {
    String genome_build
    File input_vcf
    File input_vcf_idx
    String sampleName 
    String? caller
    Array[File]? snps_vcf
    String snvIndelPASS
    File annotSVtar
    String sampleName
    String typeofAnnotation ="full"     
  }

  Int space_needed_gb = 10 + round(size(snps_vcf, "GB") + size(input_vcf, "GB")+ size(annotSVtar, "GB"))

  runtime {
    memory: "8GB"
    docker: "trinhanne/annotsv:3/1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<

    # untar the file

    tar xvzf ${annot} -C "AnnotSV"

    /opt/AnnotSV_3.1/bin/AnnotSV -bedtools /usr/bin/bedtools -outputDir "$PWD" \
    -genomeBuild ~{genome_build} -annotationDir AnnotSV -typeOfAnnotation ~{typeofAnnotation}\ 
    -SVinputFile ~{input_vcf} -snvIndelPASS ~{snvIndelPASS} \
    -outputFile ~{sampleName}.~{caller}.annotSV.tsv 
    #### -vcfFiles ~{sep="," snps_vcf}
  >>>

  output {
    File sv_variants_tsv = "~{sampleName}.~{caller}.annotSV.tsv "
  }
}
