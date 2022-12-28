## STAR fusion

version 1.0

import "https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/Terra-1.10.1/WDL/star_fusion_workflow.wdl" as starfusion

workflow STARfusion {

    input {
        String sample_id
        String runMode # Fastq or ChimJunc
        File genome_plug_n_play_tar_gz
        # Inputs required for full pipeline
        File? left_fq
        File? right_fq
        File? fastq_pair_tar_gz
        String? fusion_inspector
        Boolean examine_coding_effect
        # Inputs required for full chimeric pipeline
        File? ChimericJunction
    
        # runtime params
        String docker = "trinityctat/starfusion:latest"
        Int num_cpu = 12
        Float fastq_disk_space_multiplier = 3.25
        String memory = "50G"
        Float genome_disk_space_multiplier = 2.5
        Int preemptible = 2
        Float extra_disk_space = 10
        Boolean use_ssd = true
    }

    File ChimJun = select_first([ChimericJunction, "NULL"])

    if ( runMode == "Fastq" ){
        call starfusion.star_fusion as star_fusion {
            input: 
                left_fq = left_fq,
                right_fq = right_fq,
                fastq_pair_tar_gz = fastq_pair_tar_gz,
                genome = genome_plug_n_play_tar_gz,
                sample_id = sample_id,
                examine_coding_effect = examine_coding_effect,
                preemptible = preemptible,
                docker = docker,
                cpu = num_cpu,
                memory = memory,
                extra_disk_space = extra_disk_space,
                fastq_disk_space_multiplier = fastq_disk_space_multiplier,
                genome_disk_space_multiplier = genome_disk_space_multiplier,
                fusion_inspector = fusion_inspector,
                use_ssd = use_ssd
        }
    }

    if ( runMode == "ChimJunc" ){
        call ChimFusion {
            input: 
                ChimericJunction = ChimJun,
                genome = genome_plug_n_play_tar_gz,
                sample_id = sample_id,
                preemptible = preemptible,
                docker = docker,
                num_cpu = num_cpu,
                machine_mem_gb = memory,
                fastq_disk_space_multiplier = fastq_disk_space_multiplier,
                genome_disk_space_multiplier = genome_disk_space_multiplier,
                examine_coding_effect = examine_coding_effect
        }
    }

    output {
        File fusion_predictions = select_first([star_fusion.fusion_predictions, ChimFusion.fusions])
        File fusion_predictions_abridged = select_first([star_fusion.fusion_predictions_abridged, ChimFusion.fusionsAbridged])
        File? junction = star_fusion.junction
        File? bam = star_fusion.bam
        File? sj = star_fusion.sj
        File? coding_effect = star_fusion.coding_effect
        Array[File]? extract_fusion_reads = star_fusion.extract_fusion_reads
        File? star_log_final = star_fusion.star_log_final    
        File? fusion_inspector_validate_fusions_abridged = star_fusion.fusion_inspector_validate_fusions_abridged
        File? fusion_inspector_validate_web = star_fusion.fusion_inspector_validate_web
        File? fusion_inspector_inspect_fusions_abridged = star_fusion.fusion_inspector_inspect_fusions_abridged
        File? fusion_inspector_inspect_web = star_fusion.fusion_inspector_inspect_web
    }
}

task ChimFusion {
    input {
      File ChimericJunction
      File genome
      String docker
      String sample_id 
      Int num_cpu
      Float fastq_disk_space_multiplier
      Float genome_disk_space_multiplier
      String? machine_mem_gb
      String? preemptible
      Boolean examine_coding_effect
    }

    Int disk_space_gb = ceil(genome_disk_space_multiplier*size(genome, "GB")+fastq_disk_space_multiplier*size(ChimericJunction, "GB"))

    command <<<

    mkdir -p genome_dir

    tar xf ~{genome} -C genome_dir --strip-components 1

    # Identify the fusions from chimericjunction file

    /usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
             -J ~{ChimericJunction} \
             --output_dir ~{sample_id} --CPU ~{num_cpu} ~{true='--examine_coding_effect' false='' examine_coding_effect}

    >>>

    output {
      File fusions = "~{sample_id}/star-fusion.fusion_predictions.tsv"
      File fusionsAbridged = "~{sample_id}/star-fusion.fusion_predictions.abridged.tsv"
    }

    runtime {
      docker: docker
      memory: select_first([machine_mem_gb, 50]) + " GB"
      cpu: num_cpu
      disks: "local-disk " + disk_space_gb + " HDD"
      preemptible: select_first([preemptible, 3])
    }
}
