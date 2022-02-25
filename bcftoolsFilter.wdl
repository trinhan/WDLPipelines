version 1.0

workflow filterPASS {
	input {
		File vcfgz
		File vcftbi
		Boolean findPass
		String? SearchTerm
	}

	call runbcftools {
		input:
		vcfgz=vcfgz,
		vcftbi=vcftbi,
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
		Boolean findPass
		String? SearchTerm	
	}
		Boolean FiltTerm = if defined(SearchTerm) then true else false
		Int disk_space_gb = 2*ceil(size(vcfgz, "GB"))
		String vcf_basename = basename(vcfgz, ".vcf.gz") 

	command {


		bcftools view ${true="-f 'PASS,.'" false="" findPass} ~{"-i " + SearchTerm} ~{vcfgz} > ~{vcf_basename}.pass.vcf

		bgzip ~{vcf_basename}.pass.vcf
		tabix -p vcf ~{vcf_basename}.pass.vcf.gz
	}

	runtime {
    	docker: "trinhanne/sambcfhts:v1.13.3"
    	memory: "5 GB"
    	cpu: "1"
    	disks: "local-disk " + disk_space_gb + " HDD"
    	preemptible: 3
  }
  
  output {
    File Outvcfgz = "~{vcf_basename}.pass.vcf.gz"
    File Outvcfgztbi = "~{vcf_basename}.pass.vcf.gz.tbi"
  }
}