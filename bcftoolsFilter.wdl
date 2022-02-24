version 1.0

workflow filterPASS {
	input {
		File vcfgz
		File vcftbi
		String outname
		Boolean findPass
		String? SearchTerm
	}

	call runbcftools {
		input:
		vcfgz=vcfgz,
		vcftbi=vcftbi,
		outname=outname,
		findPass=findPass,
		SearchTerm=SearchTerm
	}

	output {
		File Outvcfgz = runbcftools.Outvcfgz
		File Outvcfgztbi = runbcftools.Outvcfgztbi
	}
}

task runbcftools{
	input {
		File vcfgz
		File vcftbi
		String outname
		Boolean findPass
		String? SearchTerm	
	}
		Boolean FiltTerm = if defined(SearchTerm) then true else false
		Int disk_space_gb = 2*ceil(size(vcfgz, "GB"))

	command {


		bcftools view ${true="-f 'PASS,.'" false="" findPass} ~{"-i " + SearchTerm} ~{vcfgz} > ~{outname}.vcf

		bgzip ~{outname}.vcf
		tabix -p vcf ~{outname}.vcf.gz
	}

	runtime {
    	docker: "trinhanne/sambcfhts:v1.13.3"
    	memory: "5 GB"
    	cpu: "1"
    	disks: "local-disk " + disk_space_gb + " HDD"
    	preemptible: 3
  }
  
  output {
    File Outvcfgz = "~{outname}.vcf.gz"
    File Outvcfgztbi = "~{outname}.vcf.gz.tbi"
  }
}