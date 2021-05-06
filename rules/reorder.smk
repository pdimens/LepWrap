rule reorder_markers:
    input:
        datacall = "2_Filtering/data_f.call.gz",
        filt_map = "map.master",
        lg_order = "5_Trim/{trimfile}.trimmed"
    output:
        "6_OrderMarkers/iterations/{trimfile}.{ITER}"
    log:
        run = "6_OrderMarkers/logs/runs/{trimfile}.{ITER}.log",
        recomb = "6_OrderMarkers/logs/recombination/{trimfile}.{ITER}.recombinations"
    message:
        """
        Reordering {input.lg_order}, iteration: {params.iteration}
        """
    params:
        dist_method = dist_method,
        eval_order="evaluateOrder=5_Trim/{trimfile}.trimmed",
        iteration = "{ITER}"
    threads: threads_per
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.eval_order} {params.dist_method} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule reorder_summary:
    input:
        expand("6_OrderMarkers/iterations/ordered.{lg}.{iter}", lg = lg_range, iter = ITER)
    output:
        like = "6_OrderMarkers/likelihoods.summary",
        recomb = "6_OrderMarkers/recombination.summary"
    message:
        """
        Iteration likelihood summary >> {output.like}
        Recombination summary >> {output.recomb}        
        """
    shell:
        """
        for LIKE in {input}; do 
            LG=$(echo $(basename $LIKE) | cut -d "." -f1,2)
            ITERUN=$(echo $LIKE | cut -d "." -f3)
            LIKELIHOOD=$(cat $LIKE | grep "likelihood = " | cut -d " " -f7)
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> {output.like}.tmp
        done
        sort {output.like}.tmp -k1,1V -k3,3nr > {output.like} && rm {output.like}.tmp
        Rscript scripts/RecombinationSummary.r 6_OrderMarkers
        """

rule best_likelihoods2:
    input:
        "6_OrderMarkers/likelihoods.summary"
    output:
        "6_OrderMarkers/best.likelihoods"
    message:
        """
        Identifying ordered maps with best likelihoods for each LG >> {output}
        """
    shell:
        """
        NUMITER=$(find 6_OrderMarkers/iterations -maxdepth 1 -name "ordered.*.*" | cut -d "." -f3 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find 6_OrderMarkers/iterations -maxdepth 1 -name "ordered.*.*" | wc -l) 
        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p {input} | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "6_OrderMarkers/iterations/$LIKELYMAP" >> {output}
        done
        """