rule trim_edge_clusters:
    input: "4_OrderMarkers/ordered.{lg_range}"
    output: "5_Trim/ordered.{lg_range}.trimmed"
    log:
        "5_Trim/logs/ordered.{lg_range}.removed",
        "5_Trim/logs/ordered.{lg_range}.trim.pdf"
    params:
        trim_threshold = trim_thresh,
        edge_length = edge_len
    message: "Removing edge clusters >{params.trim_threshold}cM apart from the other markers at the ends of {input}"
    shell:
        """
        Rscript scripts/LepWrapTrim.r {input} {params.trim_threshold} {params.edge_length}
        """

rule trim_summary:
    input: expand("5_Trim/ordered.{lg}.trimmed", lg = lg_range)
    output:
        detailed = "5_Trim/trim.details",
        summary = "5_Trim/trim.summary"
    message: "Summarizing trim logs to {output}"
    priority: 1
    shell:
        """
        echo "# this is a summary of which markers were removed from which linkage group via trimming distant edge clusters" >> {output.detailed}
        echo -e "LG\trm_marker" >> {output.detailed}
        for each in 5_Trim/logs/ordered.*.removed ; do
            BASE=$(basename $each | cut -d "." -f1,2)
            sed -e "s/^/$BASE /" $each >> {output.detailed}.tmp 
        done
        sort -V {output.detailed}.tmp > {output.detailed} && rm {output.detailed}.tmp
        echo "n_removed map" > {output.summary}.tmp
        cut -d" " -f1 {output.detailed} | uniq -c  >> {output.summary}.tmp
        column -t {output.summary}.tmp > {output.summary} && rm {output.summary}.tmp 
        """
