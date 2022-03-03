## GATK selectVariants
## script to filter variants using gatk
## useful for removing <NON_REF> arguments
##
## Inputs:
## - vcf.gz file
## - vcfgztbi: index file required if .gz file needed
## - referenceFasta: optional, only needed for some options
## - searchStrings e.g. "--exclude-non-variants --remove-unused-alternates"
##
## Outputs:
## - vcf filtered file
## - vcf index file

version 1.0

workflow runSelectVariants {
	input {
		File vcfgz
		File? vcftbi
		File? referenceFasta
		String? searchString
	}


	call SelectVariants {
		input:
		vcfgz=vcfgz,
		vcftbi=vcftbi,
		referenceFasta=referenceFasta,
		searchString=searchString
	}

	output {
		File OutvcfgzSR = SelectVariants.Outvcfgz
		File OutvcfgztbiSR = SelectVariants.Outvcfgztbi
	}

}

task SelectVariants{
	input {
		File vcfgz
		File? vcftbi
		File? referenceFasta
		String? searchString
	}

	String vcf_basename = basename(vcfgz, ".vcf.gz") 
	Int disk_space_gb = 2*ceil(size(vcfgz, "GB")+ size(referenceFasta, "GB")+1)

	command {

	cp ~{vcfgz} .
	cp ~{vcftbi} .

	gatk SelectVariants -V ~{vcf_basename}.vcf.gz ~{"-R "+ referenceFasta} ~{searchString} -O ~{vcf_basename}.SVfilt.gvcf.gz

	}

	runtime {
    	docker: "broadinstitute/gatk:4.2.4.0"
    	memory: "12 GB"
    	cpu: "1"
    	disks: "local-disk " + disk_space_gb + " HDD"
    	preemptible: 3
  }

	output {
		File Outvcfgz = "~{vcf_basename}.SVfilt.gvcf.gz"
		File Outvcfgztbi = "~{vcf_basename}.SVfilt.gvcf.gz.tbi"
	}

}