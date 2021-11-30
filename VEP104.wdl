version 1.0

## modified from
## Creator: KK
## Minor Update: 2020-07-16 BS
## scripts version: 1.0.1
##
##
# VEP on WDL
## Annotate File and produce a result format in Tab-delimited or VCF format.
##
## **VEP version**: 104
## **Cache version**: 104
##
## ** Inputs **
## species: default homo_sapiens
## input_format: default vcf
## cache: file location of cache/vep resources: default set to hg37
## plugins: plugins to use, in the form "--plugin CADD --plugin REVEL"
## OTHER_VEP_OPTS: other inputs to pass e.g. "--polyphen b"
## assembly: genome used (eg, GRCh37, GRCh38)
## refFasta is needed for HGVS annotations
## 
## Check outpouts here: [https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic]
##  
##
## VEP
## * [Download](http://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html)
## * [Run](http://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_numbers)
## * [Plugins](http://asia.ensembl.org/info/docs/tools/vep/script/vep_plugins.html)
## * [Cache](http://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
## * [Ensembl Variation - Data prediction](https://asia.ensembl.org/info/genome/variation/prediction/index.html)
##
task variant_effect_predictor {
    input {
    File inputFile
    String sample_name
    File cache
    Int? fork
    String species
    String assembly
    String input_format
    File? refFastaFai
    File? refFasta
    File? caddSnv
    File? caddSnvTbi
    File? caddIndel
    File? caddIndelTbi
    File? revelPlugin
    File? clinvar
    File? clinvarTbi
    File? gnomad
    File? gnomadIdx
    # shortcut for common flags

    # output options
    Boolean? humdiv
    # [1/0] output: is gene associated with a phenotype
    Boolean? gene_phenotype
    # overlaps with known regulatory elements
    Boolean? regulatory
    # force to interpret as phased
    Boolean? phased
    Boolean? allele_number
    # Give cDNA, CDS and protein positions as Position/Length
    Boolean? total_length
    # Adds affected exon and intron numbering to to output. Format is Number/Total
    Boolean? numbers
    # Adds names of overlapping protein domains to output
    Boolean? domains
    Boolean? no_escape
    # Don't overwrite existing CSQ entry
    Boolean? keep_csq
    Boolean? no_consequences
    Boolean? variant_class

    String? cell_type
    String? individual

    # variant prediction effects

    String? sift
    String? polyphen
    String? vcf_info_field
    String? terms

    # identifiers
    Int? shift_hgvs
    Boolean? hgvs
    Boolean? protein
    Boolean? symbol
    Boolean? ccds
    Boolean? uniprot
    # translational support level
    Boolean? tsl
    # appris isoform annotation
    Boolean? appris
    # flag for canoncial transcript
    Boolean? canonical
    # protein coding, IGR etc
    Boolean? biotype
    Boolean? xref_refseq

    # co-located variants
    Boolean? check_existing
    Boolean? check_alleles
    Boolean? check_svs
    Boolean? gmaf
    Boolean? af_1kg
    Boolean? af_esp
    Boolean? af_gnomad
    Boolean? old_maf
    Boolean? pubmed
    Int? failed

    # output format options
     Boolean? minimal

    # filtering and qc options
    Boolean? check_ref
    Boolean? coding_only
    Array[String]? chr
    Boolean? no_intergenic
    Boolean? pick
    Boolean? pick_allele
    Boolean? flag_pick
    Boolean? flag_pick_allele
    Boolean? per_gene
    Array[String]? pick_order
    Boolean? most_severe
    Boolean? summary
    Boolean? filter_common
    Boolean? check_frequency
    String? freq_pop
    Int? freq_freq
    String? freq_gt_lt
    String? freq_filter
    Boolean? allow_non_variant
    Int? buffer_size
    String? OTHER_VEP_OPTS

    ## runtime
    String docker_image = "ensemblorg/ensembl-vep:release_104.3"
    String VepDIR = "/opt/vep/src/ensembl-vep/vep"
    String VEP_DATA="/opt/vep/.vep/"
    String VEP_DATA_PLUGINS="/opt/vep/.vep/Plugins"
    
    Int cpu =select_first([fork, 4])
    Int machine_mem_gb = 4
    Int? preemptible_tries
    Int? max_retries
    }
    Int disk_space_gb = ceil(5*size(cache, "GB")+size(inputFile, "GB")+10)
    String runcadd = if (defined(caddSnv)||defined(caddIndel)) then "1" else "0"
    String runrevel = if defined(revelPlugin) then "1" else "0"
    String runclinvar = if defined(clinvar) then "1" else "0"
    String rungnomad = if defined(gnomad) then "1" else "0"

    command {
        set -e

        mkdir -v vep_data_dir
        #delete the existing directory first to make a successful link
        rm -rf ~/.vep
        ln -vs `pwd`/vep_data_dir ~/.vep

        tar xvzf ${cache} -C vep_data_dir

        plugins=""
        mkdir -v vep_data_dir/Plugins
        

        if [ ${runrevel} -eq "1" ];
        then
        unzip ${revelPlugin}
        cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
        sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
        bgzip new_tabbed_revel.tsv
        curl -O https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/104/REVEL.pm 
        mv REVEL.pm vep_data_dir/Plugins
        mv new_tabbed_revel.tsv.gz vep_data_dir/Plugins
        tabix -s 1 -b 2 -e 2 vep_data_dir/Plugins/new_tabbed_revel.tsv.gz
        plugins="$plugins --plugin REVEL,${VEP_DATA_PLUGINS}/new_tabbed_revel.tsv.gz "
        fi

        echo 'show everything'
        ls
        echo 'show everything in vep_data_dir/Plugins'
        ls vep_data_dir/Plugins
        echo 'show everything in .vep Plugins'
        ls ~/.vep/Plugins
        echo 'show everything in Plugins'
        ls /opt/vep/.vep/Plugins
        
        ##/opt/vep/src/ensembl-vep/perl INSTALL.pl -a cfp -s ${species} -y ${assembly} -g REVEL,CADD

        if [ ${runcadd} -eq "1" ];
        then
        curl -O https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/104/CADD.pm 
        mv CADD.pm vep_data_dir/Plugins
        plugins="$plugins --plugin CADD${"," + caddSnv}${"," + caddIndel} "
        fi 

        echo 'show everything in Plugins'
        ls /opt/vep/.vep/Plugins
        echo $plugins
        

        runclinvar=""

        if [ ${runclinvar} -eq "1" ];
        then
        runclinvar="$runclinar --custom ${clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN "
        fi

        if [ ${rungnomad} -eq "1" ];
        then 
        runclinvar="$runclinar --custom ${gnomad},gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH "
        fi

        ## gnomad

    
        ### Plugin
        ##    "GRCh38"    
        ## --plugin dbNSFP,/cromwell_root/data/vep/.vep/dbNSFP/4.0/dbNSFP4.0b2.gz,ALL
        ## --plugin dbscSNV,/cromwell_root/data/vep/.vep/dbscSNV/dbscSNV1.1_GRCh38.txt.gz
        ##  "GRCh37"
        ## --plugin dbNSFP,/cromwell_root/data/vep/.vep/dbNSFP/3.5/dbNSFP_hg19.gz,ALL
        ## --plugin dbscSNV,/cromwell_root/data/vep/.vep/dbscSNV/dbscSNV1.1_GRCh37.txt.gz
        ## "--plugin CADD --plugin REVEL"




        ## Running Section
        ${VepDIR} \
        --no_progress  \
        --offline --cache \
        --input_file ${inputFile} \
        --format ${input_format} \
        --species ${species} \
        --assembly ${assembly} \
        --dir_cache ${VEP_DATA} \
        --stats_file ${sample_name}.${assembly}_vep_summary.html \
        --dir_plugins ${VEP_DATA_PLUGINS} \
        ${"--fasta " + refFasta} \
        ${"--fork " + fork} \
        $plugins \
        $runclinvar \
        ${"--cell_type " + cell_type} \
        ${"--individual " + individual} \
        ${"--sift " + sift} \
        ${"--polyphen " + polyphen} \
        ${"--vcf_info_field " + vcf_info_field} \
        ${"--terms " + terms} \
        ${true="--humdiv" false="" humdiv} \
        ${true="--gene_phenotype" false="" gene_phenotype} \
        ${true="--regulatory" false="" regulatory} \
        ${true="--phased" false="" phased} \
        ${true="--allele_number" false="" allele_number} \
        ${true="--total_length" false="" total_length} \
        ${true="--numbers" false="" numbers} \
        ${true="--domains" false="" domains} \
        ${true="--no_escape" false="" no_escape} \
        ${true="--keep_csq" false="" keep_csq} \
        ${true="--no_consequences" false="" no_consequences} \
        ${true="--variant_class" false="" variant_class} \
        ${"--shift_hgvs " + shift_hgvs} \
        ${true="--hgvs" false="" hgvs} \
        ${true="--protein" false="" protein} \
        ${true="--symbol" false="" symbol} \
        ${true="--ccds" false="" ccds} \
        ${true="--uniprot" false="" uniprot} \
        ${true="--tsl" false="" tsl} \
        ${true="--appris" false="" appris} \
        ${true="--canonical" false="" canonical} \
        ${true="--biotype" false="" biotype} \
        ${true="--xref_refseq" false="" xref_refseq} \
        ${true="--check_existing" false="" check_existing} \
        ${true="--check_alleles" false="" check_alleles} \
        ${true="--check_svs" false="" check_svs} \
        ${true="--gmaf" false="" gmaf} \
        ${true="--af_1kg" false="" af_1kg} \
        ${true="--af_esp" false="" af_esp} \
        ${true="--af_gnomad" false="" af_gnomad} \
        ${true="--old_maf" false="" old_maf} \
        ${true="--pubmed" false="" pubmed} \
        ${"--failed " + failed} \
        ${true="--minimal" false="" minimal} \
        ${true="--check_ref" false="" check_ref} \
        ${true="--coding_only" false="" coding_only} \
        $CHR \
        ${true="--no_intergenic" false="" no_intergenic} \
        ${true="--pick" false="" pick} \
        ${true="--pick_allele" false="" pick_allele} \
        ${true="--flag_pick" false="" flag_pick} \
        ${true="--flag_pick_allele" false="" flag_pick_allele} \
        ${true="--per_gene" false="" per_gene} \
        $PICK_ORDER \
        ${true="--most_severe" false="" most_severe} \
        ${true="--summary" false="" summary} \
        ${true="--filter_common" false="" filter_common} \
        ${true="--check_frequency" false="" check_frequency} \
        ${"--freq_pop " + freq_pop} \
        ${"--freq_freq " + freq_freq} \
        ${"--freq_gt_lt " + freq_gt_lt} \
        ${"--freq_filter " + freq_filter} \
        ${true="--allow_non_variant" false="" allow_non_variant} \
        ${"--buffer_size " + buffer_size} \
        ${OTHER_VEP_OPTS} \
        --compress_output bgzip \
        --output_file ${sample_name}.${assembly}_vep.vcf.gz
    }
    runtime {
          docker: docker_image
          preemptible: select_first([preemptible_tries, 1])
          maxRetries: select_first([max_retries, 1])
          memory: machine_mem_gb + " GB"
          disks: "local-disk " + disk_space_gb + " HDD"
          cpu: cpu
    }
    output {
        File annotatedFile = "${sample_name}.${assembly}_vep.vcf.gz"
        File? summary_html = "${sample_name}.${assembly}_vep_summary.html"
        }

}

workflow VEP {
    input {
    File inputFile
    String sample_name
    ## Keep this string - this refers to the cache
    String species ="homo_sapiens"
    ## Keep this: default is vcf
    String input_format = "vcf"
    ## cache example "gs://kccg-cb-tools/vep104/homo_sapiens_vep_104_GRCh37.tar.gz"
    File cache 
    ###= "gs://kccg-cb-tools/vep104/CADD/GRCh37/gnomad.genomes.r2.1.1.snv.tsv.gz"
    File? caddSnv  
    ###= "gs://kccg-cb-tools/vep104/CADD/GRCh37/gnomad.genomes.r2.1.1.snv.tsv.gz.tbi"
    File? caddSnvTbi
    ## cadd example "gs://kccg-cb-tools/vep104/CADD/GRCh37/gnomad.genomes.r2.1.1.indel.tsv.gz"
    File? caddIndel 
    ## cadd example "gs://kccg-cb-tools/vep104/CADD/GRCh37/gnomad.genomes.r2.1.1.indel.tsv.gz.tbi"
    File? caddIndelTbi 
    ## example "gs://kccg-cb-tools/vep104/revel-v1.3_all_chromosomes.zip"
    File? revel
    String OTHER_VEP_OPTS
    ## eg. "GRCh37" or "GRCh38" 
    String assembly
    ## optional: use reference if determining HGVSp, HGVSc 
    File? refFasta
    File? refFastaIdx
    ## example "gs://kccg-cb-tools/vep104/ClinVar/clinvar_20211030.GRCh37.vcf.gz"
    File? clinvar
    ## example "gs://kccg-cb-tools/vep104/ClinVar/clinvar_20211030.GRCh37.vcf.gz.tbi"
    File? clinvarTbi
    ## example gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf. Probably not required, used by default in cadd
    File? gnomad
    ## example gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx
    File? gnomadIdx
    } 

    call variant_effect_predictor  {
      input:
          inputFile =inputFile,
          sample_name = sample_name,
          assembly=assembly,
          species = species ,
          input_format = input_format,
          cache=cache,
          OTHER_VEP_OPTS=OTHER_VEP_OPTS,
          caddSnv=caddSnv,
          caddSnvTbi=caddSnvTbi,
          caddIndel=caddIndel,
          caddIndelTbi=caddIndelTbi,
          revelPlugin=revel,
          refFasta = refFasta,
          refFastaFai = refFastaIdx,
          clinvar=clinvar,
          clinvarTbi=clinvarTbi,
          gnomad=gnomad,
          gnomadIdx=gnomadIdx
    }
    output {
     File VEP_annotatedFile = variant_effect_predictor.annotatedFile
     File? VEP_html = variant_effect_predictor.summary_html
    }
}
