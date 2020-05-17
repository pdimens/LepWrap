rule trim_edge_clusters:
    input:
        "ordermarkers/best.likelihoods"
    output:
        "ordermarkers/best.trim/ordered.{lg_range}.trimmed"
    log:
        "ordermarkers/logs/trimming/ordered.{lg_range}.removed",
        "ordermarkers/logs/trimming/ordered.{lg_range}.trim.pdf"
    params:
        grep_lg = "ordermarkers/iterations/ordered.{lg_range}.",
        trim_threshold = trim_thresh,
        edge_length = edge_len
    message:
        """
        Removing edge clusters >{params.trim_threshold}cM apart from the other markers in first+last 15% of {params.grep_lg}.
        """
    shell:
        """
        LG=$(grep -F {params.grep_lg} {input})
        Rscript scripts/LepMak3rQA.r $(pwd) $LG {params.trim_threshold} {params.edge_length}
        """

rule trim_summary:
    input:
        expand("ordermarkers/best.trim/ordered.{lg}.trimmed", lg = lg_range)
    output:
        "ordermarkers/trim.summary"
    message:
        "Summarizing trim logs >> {output}"
    priority: 1
    shell:
        """
        echo "# this is a summary of which markers were removed from which linkage group via trimming distant edge clusters" >> {output}
        echo -e "LG\trm_marker" >> {output}
        for each in ordermarkers/logs/trimming/ordered.*.removed ; do
            BASE=$(basename $each | cut -d "." -f1,2)
            sed -e "s/^/$BASE /" $each >> {output}.tmp 
        done
        sort -V {output}.tmp > {output} && rm {output}.tmp 
        """