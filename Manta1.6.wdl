### Manta version 1.6.0 for Paired 
### INPUTS:
### bam files (normal optional)
### reference genome data
### Intervals list optional 
### Boolean save evidence: do you want to save the evidence in a bam file?
version 1.0

task MantaSomaticSV {
    input {
    String sample_name
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index
    File ref_fasta
    File ref_fasta_index
    File? region_bed
    File? region_bed_index
    Boolean is_cram = false
    Boolean save_evidence = false
    Int disk_size
    Int cpu_num = 8 
    Int mem_gb_per_job = 1
    Int mem_size = cpu_num * mem_gb_per_job
    Int preemptible_attempts = 3
    }

    command {

        if [[ -f "${normal_bam}" ]]; then 
            normal_command_line="--normalBam ${normal_bam}"
        fi 

        python /usr/local/share/manta-1.6.0-1/bin/configManta.py --tumorBam ${tumor_bam} \
                        $normal_command_line  \
                        --referenceFasta ${ref_fasta} \
                        --runDir . --generateEvidenceBam --outputContig \
                        ${"--callRegions " + region_bed}

        python runWorkflow.py --mode local \
                         -j ${cpu_num} \
                         --memGb ${mem_size}

        # change the default names with sample prefix
        if [[ -f "${normal_bam}" ]]; then
            mv results/variants/diploidSV.vcf.gz ${sample_name}.diploidSV.vcf.gz
            mv results/variants/somaticSV.vcf.gz ${sample_name}.somaticSV.vcf.gz
        else
           touch ${sample_name}.diploidSV.vcf.gz
           touch ${sample_name}.diploidSV.vcf.gz.tbi
           mv results/variants/tumorSV.vcf.gz ${sample_name}.somaticSV.vcf.gz
        fi
        
        if [ ${save_evidence} == true ];then
        tar -zcvf ${sample_name}evidence.tar.gz results/evidence/
        fi 
        
        mv results/variants/candidateSV.vcf.gz ${sample_name}.candidateSV.vcf.gz
        mv results/variants/candidateSmallIndels.vcf.gz ${sample_name}.candidateSmallIndels.vcf.gz
       
    }
    runtime {
        docker: "quay.io/biocontainers/manta:1.6.0--h9ee0642_1"
        memory: mem_size + " GB"
        cpu: cpu_num
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_attempts
    }
    output {
        File? germline_sv_vcf = "${sample_name}.diploidSV.vcf.gz"
        File somatic_sv_vcf = "${sample_name}.somaticSV.vcf.gz"
        File candidate_sv_vcf = "${sample_name}.candidateSV.vcf.gz"
        File candidate_indel_vcf = "${sample_name}.candidateSmallIndels.vcf.gz"
        File? compressed_evidence = "${sample_name}evidence.tar.gz"
    }
}

task MantaGermline{
    input {
    String sample_name
    File? bam
    File? bam_index
    File ref_fasta
    File ref_fasta_index
    File? region_bed
    File? region_bed_index
    Boolean is_cram = false
    Boolean save_evidence = false

    Int disk_size
    Int cpu_num = 8 
    Int mem_gb_per_job = 1
    Int mem_size = cpu_num * mem_gb_per_job
    Int preemptible_attempts = 3
    }

    command {

        python /usr/local/share/manta-1.6.0-1/bin/configManta.py --bam ${bam} \
                        --referenceFasta ${ref_fasta} \
                        --runDir . --generateEvidenceBam --outputContig \
                        ${"--callRegions " + region_bed}

        python runWorkflow.py --mode local \
                         -j ${cpu_num} \
                         --memGb ${mem_size}

        # change the default names with sample prefix
        
        mv results/variants/diploidSV.vcf.gz ${sample_name}.diploidSV.vcf.gz

        if [ ${save_evidence} == true];then
        tar -zcvf ${sample_name}evidence.tar.gz results/evidence/
        fi 
    }
    runtime {
        docker: "quay.io/biocontainers/manta:1.6.0--h9ee0642_1"
        memory: mem_size + " GB"
        cpu: cpu_num
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_attempts
    }
    output {
        File germline_sv_vcf = "${sample_name}.diploidSV.vcf.gz"
        File? compressed_evidence = "${sample_name}evidence.tar.gz"
    }
}

workflow Manta {
    input {
    String sample_name
    File bam
    File bam_index
    File? normal_bam
    File? normal_bam_index
    File ref_fasta
    File ref_fasta_index
    File? region_bed
    File? region_bed_index
    Boolean save_evidence = false
    Boolean is_cram = false
    String runMode = "Tumour" ## "Tumour" or "Germline" 
    }
    
    Int disk_size=4*ceil(size(bam, "GB")+ size(normal_bam, "GB") +size(ref_fasta, "GB")) 

    if (runMode =="Tumour"){
    call MantaSomaticSV {
        input: sample_name = sample_name,
               tumor_bam = bam,
               tumor_bam_index = bam_index,
               normal_bam = normal_bam,
               normal_bam_index = normal_bam_index,
               ref_fasta = ref_fasta,
               ref_fasta_index = ref_fasta_index,
               is_cram = is_cram,
               region_bed=region_bed,
               disk_size=disk_size,
               region_bed_index=region_bed_index,
               save_evidence=save_evidence
    }
    }

    if (runMode =="Germline"){
    call MantaGermline {
        input: sample_name = sample_name,
               bam = bam,
               bam_index = bam_index,
               ref_fasta = ref_fasta,
               ref_fasta_index = ref_fasta_index,
               is_cram = is_cram,
               region_bed=region_bed,
               disk_size=disk_size,
               region_bed_index=region_bed_index,
               save_evidence=save_evidence
    }
    }

    output {
        File? germline_sv_vcf = select_first([MantaSomaticSV.germline_sv_vcf, MantaGermline.germline_sv_vcf]) 
        File? somatic_sv_vcf = MantaSomaticSV.somatic_sv_vcf
        File? candidate_sv_vcf = MantaSomaticSV.candidate_sv_vcf
        File? candidate_indel_vcf = MantaSomaticSV.candidate_indel_vcf
        File? evidence = select_first([MantaSomaticSV.compressed_evidence, MantaGermline.compressed_evidence]) 
    }
}
