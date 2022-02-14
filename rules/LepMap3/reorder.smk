rule reorder_markers:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        filt_map = "LOD.master",
        lg_order = "5_Trim/ordered.{lg_range}.trimmed"
    output: 
        lg = "6_OrderMarkers/ordered.{lg_range}",
        runlog = temp("6_OrderMarkers/logs/ordered.{lg_range}.running")
    log:
        run = "6_OrderMarkers/logs/ordered.{lg_range}.log",
        recomb = "6_OrderMarkers/recombination/ordered.{lg_range}.recombination"
    message: "Reordering linkage group {params.lg}"
    params:
        lg = "{lg_range}",
        extra = reorder_extra
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp $CONDA_PREFIX/bin/ OrderMarkers2 evaluateOrder={input.lg_order} {params.extra} map={input.filt_map} data=- numThreads={threads} &> {output.runlog}
        sed -n '/\*\*\* LG \=/,$p' {output.runlog} > {output.lg}
        grep "recombin" {output.runlog} > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {output.runlog} > {log.run}
        """

rule reorder_summary:
    input: expand("6_OrderMarkers/ordered.{lg}", lg = lg_range)
    output: "6_OrderMarkers/recombination/recombination.summary"
    message: "Recombination summary of reordering: {output}"
    shell:
        """
        RecombinationSummary.r 6_OrderMarkers/recombination > {output}
        """
