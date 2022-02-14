rule calculate_distances:
    input:
        data_call = "2_Filtering/data.filtered.lepmap3.gz",
        lg = "6_OrderMarkers/ordered.{lg_range}"
    output:
        distance = "7_Distances/ordered.{lg_range}.distances",
        sex_averaged = "7_DistancesSexAverage/ordered.{lg_range}.sexavg",
        sex_averagedtmp = temp("7_DistancesSexAverage/logs/.ordered.{lg_range}.sexavg"),
        intervals = "7_Intervals/ordered.{lg_range}.intervals"
    message: "Calculating marker distances and intervals for linkage group: {params.lg}"
    log:
        sex_averaged = "7_DistancesSexAverage/logs/ordered.{lg_range}.sexavg.log",
        intervals = "7_Intervals/logs/ordered.{lg_range}.intervals.log"
    params:
        dist_method = dist_method,
        lg = "{lg_range}"
    threads: 2
    shell:
        """
        cp {input.lg} {output.distance}
        
        zcat {input.data_call} | java -cp $CONDA_PREFIX/bin/lepmap3 OrderMarkers2 evaluateOrder={input.lg} data=- {params.dist_method} numThreads={threads} improveOrder=0 sexAveraged=1 &> {output.sex_averagedtmp}
        sed -i -e 's/LG \= 0/LG \= {params.lg}/g' {output.sex_averagedtmp}
        sed -n '/\*\*\* LG \=/,$p' {output.sex_averagedtmp} > {output.sex_averaged} 
        awk '/#java/{{flag=1}} flag; /*** LG =/{{flag=0}}' {output.sex_averagedtmp} > {log.sex_averaged}

        zcat {input.data_call} | java -cp $CONDA_PREFIX/bin/lepmap3 OrderMarkers2 evaluateOrder={input.lg} data=- {params.dist_method} numThreads={threads} calculateIntervals={output.intervals} > {log.intervals} 2>&1
        """