## Workflow: Header, Map Rate, Meandepth
## Inputs:
##   (I1) BAM File (.bam)
##   (I2) BAM Index (.bai)

## Outputs:
##   (O1) File containing the header of the BAM file
##   (O2) File indicating the mapping rate
##   (O3) File indicating the average depth

version 1.0
workflow HeaderMaprateDepth{
	input{
		File input_bam
		File index_bam
		String filename
		Int preemptible = 3
		Int cpu = 1
		Int memoryGB = 8
	}

	call mapCaller{
		input:
			input_bam=input_bam,
			index_bam=index_bam,
			filename=filename,
			memoryGB=memoryGB,
			preemptible=preemptible,
			cpu=cpu
	}

	call depthCaller{
		input:
			input_bam=input_bam,
			index_bam=index_bam,
			filename=filename,
			memoryGB=memoryGB,
			preemptible=preemptible,
			cpu=cpu
	}
 
	output {
		File maprate = mapCaller.maprate
		File meandepth = depthCaller.meandepth
	}
}


# This task calls the header from a BAM file.
task headCaller{
	input{
		File input_bam
		File index_bam
		String filename
		Int preemptible
		Int cpu
		Int memoryGB
	}

	Int diskSpace = 3*ceil(size(input_bam, "GB") + size(index_bam, "GB"))

	command {
		samtools view -H ${input_bam} > ${filename}_header.txt
	}

	output {
		File header = "${filename}_header.txt"
	}

	runtime {
		docker: "staphb/samtools:1.14"
    	preemptible: preemptible
    	disks: "local-disk ${diskSpace} HDD"
    	cpu: cpu
    	memory: memoryGB
	}
}
# This task calls the mapping rate from a BAM file by highlighting the "%" chr.
task mapCaller {
	input {
		File input_bam
		File index_bam
		String filename
		Int preemptible
		Int cpu
		Int memoryGB
	}
	
	Int diskSpace = 3*ceil(size(input_bam, "GB") + size(index_bam, "GB"))

	command {
		samtools flagstat ${input_bam} | grep "%" > ${filename}_maprate.txt
	}

	output {
		File maprate = "${filename}_maprate.txt"
	}

	runtime {
		docker: "staphb/samtools:1.14"
    	preemptible: preemptible
    	disks: "local-disk ${diskSpace} HDD"
    	cpu: cpu
    	memory: memoryGB
	}
}
# This task retrieves the mean depth from a BAM file by means of samtools:coverage.
task depthCaller {
	input {
		File input_bam
		File index_bam
		String filename
		Int preemptible
		Int cpu
		Int memoryGB
	}
	
	Int diskSpace = 3*ceil(size(input_bam, "GB") + size(index_bam, "GB"))

	command {
		samtools coverage -m ${input_bam} > ${filename}_meandepth.txt
	}

	output {
		File meandepth = "${filename}_meandepth.txt"
	}

	runtime {
		docker: "staphb/samtools:1.14"
    	preemptible: preemptible
    	disks: "local-disk ${diskSpace} HDD"
    	cpu: cpu
    	memory: memoryGB
	}
}
