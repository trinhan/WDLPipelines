version 1.0

workflow downsampling {

    input {
	File file_bam
	File file_bai
    }
    
    call downsampling {
	input:
	file_bam = file_bam,
	file_bai = file_bai
    }

    output {
	File downsampling_output = downsampling.file_out
    }

}	 

task downsampling {

    input {
	File file_bam
	File file_bai
	Int interval = 1000000
	String? tag
	String roi
	String docker_image
	Int memory_size = 32
    }

    String base_name = basename(file_bam, ".bam")
    String file_out = base_name + ".txt"
    
    command {
	Rscript /opt/downsampling_1.1.R --args \
	--bam=${file_bam} \
	--out=${file_out} \
	--int=${interval} \
	${"--tag=" + tag} \
	${"--roi=" + roi}
    }

    output {
 	File file_out = "${file_out}"
    }
    
    Int disk_size = 2 * ceil(size(file_bam, "GB"))

    runtime {
	docker: docker_image
	memory: "${memory_size} GB"
	disks: "local-disk ${disk_size} HDD"
    }

}
