rule order_markers:
    input:
        datacall = "2_Filtering/data_f.call.gz",
        filt_map = "map.master"
    output:
        "4_OrderMarkers/iterations/ordered.{lg_range}.{ITER}"
    log:
        run = "4_OrderMarkers/logs/runs/ordered.{lg_range}.{ITER}.log",
        recomb = "4_OrderMarkers/logs/recombination/ordered.{lg_range}.{ITER}.recombinations"
    message: "Ordering the markers with {params.dist_method} on linkage group: {params.chrom}, iteration: {params.iteration}"
    params:
        dist_method = dist_method,
        chrom = "{lg_range}",
        iteration = "{ITER}"
    threads: threads_per
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.dist_method} chromosome={params.chrom} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule order_summary:
    input:
        expand("4_OrderMarkers/iterations/ordered.{lg}.{iter}", lg = lg_range, iter= ITER)
    output:
        like = "4_OrderMarkers/likelihoods.summary",
        recomb = "4_OrderMarkers/recombination.summary"
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
        Rscript scripts/RecombinationSummary.r 4_OrderMarkers
        """

rule best_likelihoods:
    input:
        "4_OrderMarkers/likelihoods.summary"
    output:
        "4_OrderMarkers/best.likelihoods"
    message:
        """
        Identifying orders with best likelihoods for each LG >> {output}
        """
    shell:
        """
        NUMITER=$(find 4_OrderMarkers/iterations -maxdepth 1 -name "ordered.*.*" | cut -d "." -f3 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find 4_OrderMarkers/iterations -maxdepth 1 -name "ordered.*.*" | wc -l) 
        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p {input} | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "4_OrderMarkers/iterations/$LIKELYMAP" >> {output}
        done
        """