rule reorder:
    input:
        datacall = "data_f.call.gz",
        filt_map = "map.master",
        lg_order = "ordermarkers/best.trim/{trimfile}.trimmed"
    output:
        "reordermarkers/iterations/{trimfile}.{ITER}"
    log:
        run = "reordermarkers/logs/runs/{trimfile}.{ITER}.log",
        recomb = "reordermarkers/logs/recombination/{trimfile}.{ITER}.recombinations"
    message:
        """
        Reordering {input.lg_order}, iteration: {params.iteration}
        """
    params:
        dist_method = "useKosambi=1",
        eval_order="evaluateOrder=ordermarkers/best.trim/{trimfile}.trimmed",
        iteration = "{ITER}"
    threads: config["threads_per"]
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.eval_order} {params.dist_method} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule reorder_summary:
    input:
        expand("reordermarkers/iterations/ordered.{LG}.{iter}", LG = lg_range, iter = ITER)
    output:
        like = "reordermarkers/likelihoods.summary",
        recomb = "reordermarkers/recombination.summary"
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
        Rscript scripts/RecombinationSummary.r reordermarkers
        """

rule find_bestlikelihoods2:
    input:
        "reordermarkers/likelihoods.summary"
    output:
        "reordermarkers/best.likelihoods"
    message:
        """
        Identifying ordered maps with best likelihoods for each LG >> {output}
        """
    shell:
        """
        NUMITER=$(find reordermarkers/iterations -maxdepth 1 -name "ordered.*.*" | cut -d "." -f3 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find reordermarkers/iterations -maxdepth 1 -name "ordered.*.*" | wc -l) 
        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p {input} | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "reordermarkers/iterations/$LIKELYMAP" >> {output}
        done
        """