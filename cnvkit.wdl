version 1.0

workflow cnvkit {
    input {
    File input_bam
    File input_bai
    File? reference_cnn
    File target_bed
    File? antitarget_bed
    }

    call cnvkit_coverage {
    	input:
            input_bam = input_bam,
            input_bai = input_bai,
            target_bed=target_bed,
            antitarget_bed=antitarget_bed
    }

    if (defined(reference_cnn)) {  
    call cnvkit_analysis {
        input:
            input_target_cnn = cnvkit_coverage.output_target_cnn,
            input_antitarget_cnn = cnvkit_coverage.output_antitarget_cnn,
            reference_cnn = reference_cnn
    }
    }

    output {
    File cnvkit_target_cnn = cnvkit_coverage.output_target_cnn
    File cnvkit_antitarget_cnn = cnvkit_coverage.output_antitarget_cnn
    File? cnvkit_cnr = cnvkit_analysis.output_cnr
    File? cnvkit_cns = cnvkit_analysis.output_cns
    File? cnvkit_seg = cnvkit_analysis.output_seg
    File? cnvkit_scatter = cnvkit_analysis.output_scatter
    File? cnvkit_diagram = cnvkit_analysis.output_diagram        
    File? cnvkit_breaks = cnvkit_analysis.output_breaks        
    File? cnvkit_genemetrics = cnvkit_analysis.output_genemetrics        
    File? cnvkit_metrics = cnvkit_analysis.output_metrics
    File? antitarget_bedout = cnvkit_coverage.antitarget_bedOut 
    }
}     

task cnvkit_analysis {
    input {
    File input_target_cnn
    File input_antitarget_cnn
    File? reference_cnn
    String docker_image = "etal/cnvkit:0.9.8"
    Int memory_size = 16
    Int threads = 4
    Int preemptible_tries = 3
    Int? disk_size = 100
    }


    String base_name = basename(input_target_cnn, ".targetcoverage.cnn")
    String output_cnr = base_name + ".cnr"
    String output_cns = base_name + ".cns"
    String output_seg = base_name + ".seg"
    String output_scatter = base_name + ".scatter.pdf"
    String output_diagram = base_name + ".diagram.pdf"
    String output_breaks = base_name + ".breaks.txt"
    String output_genemetrics = base_name + ".genemetrics.txt"
    String output_metrics = base_name + ".metrics.txt"

    command <<<
        cnvkit.py fix ${input_target_cnn} ${input_antitarget_cnn} ${reference_cnn} -o ${output_cnr}
        cnvkit.py segment ${output_cnr} -o ${output_cns} -p ${threads}
        cnvkit.py export seg ${output_cns} -o ${output_seg}
        cnvkit.py scatter ${output_cnr} -s ${output_cns} -o ${output_scatter}
        cnvkit.py diagram ${output_cnr} -s ${output_cns} -o ${output_diagram} -x male
        cnvkit.py breaks ${output_cnr} ${output_cns} -o ${output_breaks}
        cnvkit.py genemetrics ${output_cnr} -s ${output_cns} -o ${output_genemetrics} -x male
        cnvkit.py metrics ${output_cnr} -s ${output_cns} -o ${output_metrics}
    
    >>>

    output {
    File output_cnr = "${output_cnr}"
    File output_cns = "${output_cns}"
    File output_seg = "${output_seg}"
    File output_scatter = "${output_scatter}" 
    File output_diagram = "${output_diagram}"
    File output_breaks = "${output_breaks}" 
    File output_genemetrics = "${output_genemetrics}" 
    File output_metrics = "${output_metrics}" 
    }

    runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "${memory_size} GB"
    disks: "local-disk ${disk_size} HDD"
    }
}


task cnvkit_coverage {
    input {
        File input_bam
        File input_bai
        File target_bed
        File? antitarget_bed
        String docker_image = "etal/cnvkit:0.9.8"
        Int memory_size = 16
        Int threads = 4
        Int preemptible_tries = 3
    }
        String base_name = basename(input_bam, ".bam")
        String output_target_cnn = base_name + ".targetcoverage.cnn"
        String output_antitarget_cnn = base_name + ".antitargetcoverage.cnn"
        Boolean createAntitarget = if defined(antitarget_bed) then 0 else 1
        Int disk_size = 2 * ceil(size(input_bam, "GB"))

    command <<<
    cnvkit.py coverage ${input_bam} ${target_bed} -o ${output_target_cnn} -p ${threads}

    if [ ~{createAntitarget} -eq 1 ];
    then 
        cnvkit.py antitarget ~{target_bed} -g data/access-5kb-mappable.hg19.bed -o my_antitargets.bed
        cnvkit.py coverage ${input_bam} my_antitargets.bed -o ${output_antitarget_cnn} -p ${threads}
    else 
        cnvkit.py coverage ${input_bam} ~{antitarget_bed} -o ${output_antitarget_cnn} -p ${threads}
    fi
    >>>

    output {
        File output_target_cnn = "${output_target_cnn}"
        File output_antitarget_cnn = "${output_antitarget_cnn}"
        File? antitarget_bedOut = "my_antitargets.bed"
    }
    runtime {
        preemptible: preemptible_tries
        docker: docker_image
        memory: "${memory_size} GB"
        disks: "local-disk ${disk_size} HDD"
    }
}