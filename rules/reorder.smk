rule reorder_markers:
    input:
        datacall = "2_Filtering/data_f.call.gz",
        filt_map = "map.master",
        lg_order = "5_Trim/{trimfile}.trimmed"
    output:
        "6_OrderMarkers/{trimfile}"
    log:
        run = "6_OrderMarkers/logs/{trimfile}.{log",
        recomb = "6_OrderMarkers/recombination/{trimfile}.recombination"
    message: "Reordering {input.lg_order} with {params.iterations} iterations"
    params:
        dist_method = dist_method,
        eval_order="evaluateOrder=5_Trim/{trimfile}.trimmed",
        iterations = ITER
    threads: threads_per
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.eval_order} {params.dist_method} numMergeIterations={params.iterations} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule reorder_summary:
    input:
        expand("6_OrderMarkers/iterations/ordered.{lg}.{iter}", lg = lg_range, iter = ITER)
    output:
        recomb = "6_OrderMarkers/recombination.summary"
    message: "Reordering recombination summary: {output.recomb}"
    shell:
        """
        Rscript scripts/RecombinationSummary.r 6_OrderMarkers/recombination
        """