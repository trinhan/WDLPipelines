version 1.0

task PassJudgement {
	input {
		File? blat_reject
		File? abra2_calls
		File vcf
    	File? gnomad
    	File? m2_pon
    	String pairName
    	String refGenome
    	String runMode
    	Boolean runtest

    	# RUNTIME INPUT PARAMS
    	String? preemptible = "2"
    	String? diskGB_boot = "10"
    	String? diskGB_buffer ="5"
    	String? memoryGB ="5"
    	String? cpu ="1"
	}

	# DEFAULT VALUES
       
    Int diskGB = ceil(size(m2_pon, "G")+size(gnomad, "G")+size(blat_reject, "G"))*3 +  diskGB_buffer
    String rTest = if (runtest == true ) then "1" else "0"

	command {
		# Use an R script for this?
		# change blat to vcf file
		# bcftools isec for PoN

		if [ ${rTest} -eq "1" ];
		then
		Rscript /app/MergeVCFRecall.R --vcffile ${vcf} --outputstr ${pairName} --genome ${refGenome} --runMode ${runMode} ${"--blat " + blat_reject} ${"--pon " + m2_pon} ${"--gnomad " + gnomad} ${"--abra2 " + abra2_calls} 
		fi
	}

	runtime {
        docker         : "trinhanne/ngs_renv:0.1"       
        bootDiskSizeGb : diskGB_boot 
        preemptible    : preemptible
        cpu            : cpu
        disks          : "local-disk ~{diskGB} HDD"
        memory         : memoryGB + "GB"
    }

    output {
    File pass_judgement="${pairName}.judgement.vcf.gz"
    }
}

workflow testJudgement {
	input {
		File? blat_reject
		File? abra2_calls
		File vcf
    	File? gnomad
    	File? m2_pon
    	String pairName
    	String refGenome
    	String runMode
    	Boolean runtest
	}

	call PassJudgement {
		input:
			blat_reject=blat_reject,
			abra2_calls=abra2_calls,
			vcf=vcf,
    		gnomad=gnomad,
    		m2_pon=m2_pon,
    		pairName=pairName,
    		refGenome=refGenome,
    		runMode=runMode,
    		runtest=runtest
	}
}