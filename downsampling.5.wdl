version 1.0

workflow downsamplingWF {
    input {
	File file_bam
	File file_bai
	Int interval = 1000000
	String? tag
	String? roi
	Int memory_size = 32
    }
    
    call downsamplingTask {
	input:
	file_bam = file_bam,
	file_bai = file_bai,
	interval=interval,
	tag=tag,
	roi=roi,
	memory_size=memory_size
    }

    output {
	File downsampling_txtfile = downsamplingTask.txtout
    }

}	 

task downsamplingTask {
    input {
	File file_bam
	File file_bai
	Int interval = 1000000
	String? tag
	String? roi
	Int memory_size = 32
    }

    String base_name = basename(file_bam, ".bam")
    String file_out = base_name + ".txt"
    
    command {
	Rscript /opt/downsampling_1.2.R --args \
	--bam=${file_bam} \
	--out=${file_out} \
	--int=${interval} \
	${"--tag=" + tag} \
	${"--roi=" + roi}
    }

    output {
 	File txtout = "${file_out}"
    }
    
    Int disk_size = 2 * ceil(size(file_bam, "GB"))

    runtime {
	docker: "trinhanne/downsamplebam:v1"
	memory: "${memory_size} GB"
	disks: "local-disk ${disk_size} HDD"
    }

}
