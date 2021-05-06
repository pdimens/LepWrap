rule trim_edge_clusters:
    input:
        "4_OrderMarkers/best.likelihoods"
    output:
        "5_Trim/ordered.{lg_range}.trimmed"
    log:
        "5_Trim/logs/ordered.{lg_range}.removed",
        "5_Trim/logs/ordered.{lg_range}.trim.pdf"
    params:
        grep_lg = "4_OrderMarkers/iterations/ordered.{lg_range}.",
        trim_threshold = trim_thresh,
        edge_length = edge_len
    message:
        """
        Removing edge clusters >{params.trim_threshold}cM apart from the other markers in first+last 15% of {params.grep_lg}.
        """
    shell:
        """
        LG=$(grep -F {params.grep_lg} {input})
        Rscript scripts/LepWrapQA.r $(pwd) $LG {params.trim_threshold} {params.edge_length}
        """

rule trim_summary:
    input:
        expand("5_Trim/ordered.{lg}.trimmed", lg = lg_range)
    output:
        "5_Trim/trim.summary"
    message:
        "Summarizing trim logs >> {output}"
    priority: 1
    shell:
        """
        echo "# this is a summary of which markers were removed from which linkage group via trimming distant edge clusters" >> {output}
        echo -e "LG\trm_marker" >> {output}
        for each in 5_Trim/logs/ordered.*.removed ; do
            BASE=$(basename $each | cut -d "." -f1,2)
            sed -e "s/^/$BASE /" $each >> {output}.tmp 
        done
        sort -V {output}.tmp > {output} && rm {output}.tmp 
        """
