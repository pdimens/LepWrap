rule trim_edge_clusters:
    input: "4_OrderMarkers/ordered.{lg_range}"
    output: "5_Trim/ordered.{lg_range}.trimmed"
    log:
        "5_Trim/logs/ordered.{lg_range}.removed",
        "5_Trim/plots/ordered.{lg_range}.trim.pdf"
    params:
        trim_threshold = trim_thresh,
        edge_length = edge_len
    message: "Removing edge clusters >{params.trim_threshold}cM apart from the other markers at the ends of {input}"
    shell:
        """
        Rscript scripts/LepWrapTrim.r {input} {params.trim_threshold} {params.edge_length}
        """

rule trim_summary:
    input: 
        lg = expand("5_Trim/ordered.{lg}.trimmed", lg = lg_range),
        plots = expand("5_Trim/plots/ordered.{lg}.trim.pdf", lg = lg_range)
    output:
        detailed = "5_Trim/trim.details",
        summary = "5_Trim/trim.summary",
        summarypdf = "5_Trim/trim.summary.pdf",
        summarysvg = "5_Trim/trim.summary.svg",
        mergeplots = "5_Trim/plots/all.trimplots.pdf"
    message: "Summarizing trimming results"
    params:
        lg = lg_count
    priority: 1
    shell:
        """
        for each in 5_Trim/logs/ordered.*.removed ; do
            BASE=$(basename $each | cut -d "." -f1,2)
            sed -e "s/^/$BASE /" $each >> {output.detailed}.tmp
        done | sort -V > {output.detailed}
        #sort -V {output.detailed}.tmp >> {output.detailed} && rm {output.detailed}.tmp
        scripts/TrimCounts.r {output.detailed} {params.lg} > {output.summary}
        scripts/TrimSummaryPlot.r {output.summary}
        echo "Merging QC plots for all linkage groups"
        convert -density 300 {input.plots} {output.mergeplots}
        """
