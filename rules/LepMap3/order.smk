rule order_markers:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        filt_map = "LOD.master"
    output: 
        lg = "4_OrderMarkers/ordered.{lg_range}",
        runlog = temp("4_OrderMarkers/logs/ordered.{lg_range}.running")
    log:
        run = "4_OrderMarkers/logs/ordered.{lg_range}.log",
        recomb = "4_OrderMarkers/recombination/ordered.{lg_range}.recombinations"
    message: "Ordering linkage group {params.chrom}"
    params:
        chrom = "{lg_range}",
        extra = order_extra
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp $CONDA_PREFIX/bin/lepmap3 OrderMarkers2 chromosome={params.chrom} map={input.filt_map} {params.extra} data=- numThreads={threads} &> {output.runlog}
        sed -n '/\*\*\* LG \=/,$p' {output.runlog} > {output.lg}
        grep "recombin" {output.runlog} > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {output.runlog} > {log.run}
        """

rule recomb_summary:
    input: expand("4_OrderMarkers/ordered.{lg}", lg = lg_range)
    output: "4_OrderMarkers/recombination/recombination.summary"
    message: "Recombination summary: {output}"
    shell: "RecombinationSummary.r 4_OrderMarkers/recombination > {output}"
