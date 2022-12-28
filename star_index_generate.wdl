## WF to generate star index (if needed)
## Requirements:
## - genomeFasta (from terra)
## - name of output index file (e.g. STAR_hg38_index). Will be compressed as a tar.gz file that needs to be unzipped in STAR pipeline
## - gencode annotation file (https://www.gencodegenes.org/human/ - using first option here, saved in kccg-cb-tools)
## Resource requirements: (index works for hg38 with the following parameters)
## - 100GB of disk space
## - 10 threads
## - 40GB memory

version 1.0

workflow star_index {
    call star_generate_index
}


task star_generate_index {
    input {
        File genomeFile
        Int nthreads = 10
        String OutputFile
        String? mem
        String? disk_space
        Int? preemptible
        File annotation
        String docker

    }

        String disk_space = select_first([disk_space, "100"])
        String memory=select_first([mem, "40"])

    command {
        set -euo pipefail
        
        mkdir star_index
        
        zcat ~{annotation} > annotation.gtf

        STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ~{genomeFile} --runThreadN ~{nthreads} --sjdbGTFfile annotation.gtf

        tar -czvf ~{OutputFile}.tar.gz star_index/
    }

    output {
        File index_file = "~{OutputFile}.tar.gz"
    }

     runtime {
        docker: select_first([docker, "broadinstitute/gtex_rnaseq:V10"])
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "~{nthreads}"
        preemptible: select_first([preemptible, 3])
    }
}