rule order_markers:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        filt_map = "LOD.master"
    output: "4_OrderMarkers/ordered.{lg_range}"
    log:
        run = "4_OrderMarkers/logs/ordered.{lg_range}.log",
        recomb = "4_OrderMarkers/recombination/ordered.{lg_range}.recombinations"
    message: "Ordering the markers with {params.dist_method} on linkage group {params.chrom}"
    params:
        chrom = "{lg_range}",
        iterations = ITER,
        phase = phasenum,
        extra = order_extra
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp software/LepMap3 OrderMarkers2 map={input.filt_map} {params.extra} data=- numThreads={threads} numMergeIterations={params.iterations} chromosome={params.chrom} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule recomb_summary:
    input:
        expand("4_OrderMarkers/ordered.{lg}", lg = lg_range)
    output:
        recomb = "4_OrderMarkers/recombination/recombination.summary"
    message: "Recombination summary: {output.recomb}"
    shell:
        """
        Rscript scripts/RecombinationSummary.r 4_OrderMarkers/recombination
        """