rule order_markers:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        filt_map = "LOD.master"
    output: "4_OrderMarkers/ordered.{lg_range}"
    log:
        runlog = temp("4_OrderMarkers/logs/ordered.{lg_range}.running"),
        run = "4_OrderMarkers/logs/ordered.{lg_range}.log",
        recomb = "4_OrderMarkers/recombination/ordered.{lg_range}.recombinations"
    message: "Ordering linkage group {params.lg} with {params.iterations} iterations"
    params:
        chrom = "{lg_range}",
        iterations = ITER,
        extra = order_extra
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp software/LepMap3 OrderMarkers2 map={input.filt_map} {params.extra} data=- numThreads={threads} numMergeIterations={params.iterations} chromosome={params.chrom} &> {log.runlog}
        sed -n '/\*\*\* LG \=/,$p' {log.runlog} > {output}
        grep "recombin" {log.runlog} > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.runlog} > {log.run}
        """

rule recomb_summary:
    input: expand("4_OrderMarkers/ordered.{lg}", lg = lg_range)
    output: "4_OrderMarkers/recombination/recombination.summary"
    message: "Recombination summary: {output}"
    shell:
        """
        Rscript scripts/RecombinationSummary.r 4_OrderMarkers/recombination > {output}
        """