rule reorder_markers:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        filt_map = "LOD.master",
        lg_order = "5_Trim/ordered.{lg_range}.trimmed"
    output:
        "6_OrderMarkers/ordered.{lg_range}"
    log:
        run = "6_OrderMarkers/logs/ordered.{lg_range}.log",
        recomb = "6_OrderMarkers/recombination/ordered.{lg_range}.recombination"
    message: "Reordering linkage group {params.lg} with {params.iterations} iterations"
    params:
        lg = "{lg_range}",
        iterations = ITER,
        extra = reorder_extra
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp software/LepMap3 OrderMarkers2 {params.extra} map={input.filt_map} data=- numThreads={threads} evaluateOrder={input.lg_order} numMergeIterations={params.iterations} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule reorder_summary:
    input:
        expand("6_OrderMarkers/ordered.{lg}", lg = lg_range)
    output:
        recomb = "6_OrderMarkers/recombination/recombination.summary"
    message: "Recombination summary of reordering: {output.recomb}"
    shell:
        """
        Rscript scripts/RecombinationSummary.r 6_OrderMarkers/recombination
        """