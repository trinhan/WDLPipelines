version 1.0
## workflow for titan
##workflow titanchunk {
##    input {
##    File tumor_hets
##    File denoised_CR
##    String pairName
##    }
##
##     call runtitan {
##        input:
##        tumor_hets = tumor_hets,
##        denoised_CR = denoised_CR,
##        pairName = pairName
##   }
##
##    output {
##        File titan_cluster_figures=runtitan.cluster_figures
##        File titan_opt_params=runtitan.opt_params
##    }
##}

task runtitan {
    input {
    File tumor_hets
    File denoised_CR
    String pairName
    String mem =3
    Int numCores =1
    Int preemptible =3
    String genomeBuild = "hg19"
    }

    command <<<
        set -e

        cat ~{tumor_hets} | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4}' | tr " " "\t" > tmp.hets.tsv
        cat ~{denoised_CR} | awk -F\t 'BEGIN{print"chr\tstart\tend\tlog2_TNratio_corrected"} /^[^@C]/ { print $0 }' > tmp.cr.tsv
        
        ## Run this for 1:3 num clusters
        
        numClusters=3
        
        ls /TitanCNA
        
        ## run TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
        echo "Maximum number of clusters: $numClusters"
        for ploidy in $(seq 2 4)
              do
            echo "Running TITAN for $i clusters."
            outDir="run_ploidy$ploidy"
            mkdir $outDir
            for numClust in $(seq 1 $numClusters)
                do
                  echo "Running for ploidy=$ploidy";
                  Rscript /TitanCNA/scripts/titanCNA.R --id ~{pairName} --hetFile tmp.hets.tsv --cnFile tmp.cr.tsv --numClusters $numClust --numCores ~{numCores} --normal_0 0.5 --ploidy_0 $ploidy --genomeBuild ~{genomeBuild} \
                  --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir $outDir --libdir /TitanCNA
                done
            echo "Completed job for $numClust clusters."
              done

        Rscript /TitanCNA/scripts/selectSolution.R --ploidyRun2=run_ploidy2 --ploidyRun3=run_ploidy3 --ploidyRun4=run_ploidy4 --threshold=0.05 --outFile optimalClusters.txt
      
        zip -r ~{pairName}_titan_figures.zip run_ploidy2 run_ploidy3 run_ploidy4

    >>>

    runtime {
    	docker: "trinhanne/titancna:1.28.0"
        cpu: numCores
        preemptible: preemptible
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
    }

    output {
        File cluster_figures = "${pairName}_titan_figures.zip"
        File opt_params = "optimalClusters.txt"
    }
}

