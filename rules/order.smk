rule order_markers:
    input:
        datacall = "2_Filtering/data_f.call.gz",
        filt_map = "map.master"
    output:
        "4_OrderMarkers/ordered.{lg_range}"
    log:
        run = "4_OrderMarkers/logs/ordered.{lg_range}.log",
        recomb = "4_OrderMarkers/recombination/ordered.{lg_range}.recombinations"
    message: "Ordering the markers with {params.dist_method} on linkage group {params.chrom}"
    params:
        dist_method = dist_method,
        chrom = "{lg_range}",
        iterations = ITER,
        phase = phasenum
    threads: threads_per
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} numMergeIterations={params.iterations} {params.dist_method} chromosome={params.chrom} phasingIterations={params.phase} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule recomb_summary:
    input:
        expand("4_OrderMarkers/ordered.{lg}", lg = lg_range)
    output:
        recomb = "4_OrderMarkers/recombination.summary"
    message: "Recombination summary: {output.recomb}"
    shell:
        """
        Rscript scripts/RecombinationSummary.r 4_OrderMarkers/recombination
        """