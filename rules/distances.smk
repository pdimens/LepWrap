rule calculate_distances:
    input:
        data_call = "data_f.call.gz",
        lg = "reordermarkers/best.likelihoods"
    output:
        distance = "distances/ordered.{lg_range}.distances",
        sex_averaged = "distances_sexaverage/ordered.{lg_range}.sexavg",
        intervals = "intervals/ordered.{lg_range}.intervals"
    message:
        """
        Calculating sex-averaged marker distances and intervals for linkage group: {params.lg}
        """
    log:
        sex_averaged = "distances_sexaverage/logs/ordered.{lg_range}.sexavg.log",
        intervals = "intervals/logs/ordered.{lg_range}.int.log"
    params:
        dist_method = "{dist_method}",
        lg = "{lg_range}"
    threads: "{threads_per}"
    shell:
        """
        LG=$(grep -F "reordermarkers/iterations/ordered.{params.lg}." {input.lg})
        cp $LG {output.distance}
        
        zcat {input.data_call} | java -cp LM3 OrderMarkers2 data=- evaluateOrder=$LG {params.dist_method} numThreads={threads} improveOrder=0 sexAveraged=1 &> {log.sex_averaged}.tmp
        sed -i -e 's/LG \= 0/LG \= {params.lg}/g' {log.sex_averaged}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.sex_averaged}.tmp > {output.sex_averaged} 
        awk '/#java/{{flag=1}} flag; /*** LG =/{{flag=0}}' {log.sex_averaged}.tmp > {log.sex_averaged} && rm {log.sex_averaged}.tmp

        zcat {input.data_call} | java -cp LM3 OrderMarkers2 data=- evaluateOrder=$LG {params.dist_method} numThreads={threads} calculateIntervals={output.intervals} > {log.intervals} 2>&1
        """