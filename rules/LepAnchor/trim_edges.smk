rule generate_updated_intervals:
  input: "13_MareyMapsUntrimmed/data.marey.gz"
  output: "14_NewIntervals/LA.intervals.{lg_range}"
  message: "Splitting out LG {params.chrom} from {input}"
  params:
    chrom = "{lg_range}"
  shell:
    """
    zgrep "LG{params.chrom}\s" {input} > {output}
    """


rule trim_newintervals:
  input: "14_NewIntervals/LA.intervals.{lg_range}"
  output: 
    outfile = "15_Trim/LA.intervals.{lg_range}.trimmed",
    plot = "15_Trim/plots/LA.intervals.{lg_range}.trim.pdf"
  message: "Trimming edge clusters for {input}"
  params:
    edge = edgelen,
    dist = trimdist
  shell: "LepWrapTrim.r {input} {params.dist} {params.edge} 15_Trim"


rule merge_trimplots:
  input: expand("15_Trim/plots/LA.intervals.{lg}.trim.pdf", lg = lg_range)
  output: "15_Trim/LA.trim.summary.pdf"
  message: "Merging trimming plots into {output}"
  shell: "convert -density 200 {input} {output}"


rule merge_trimmedintervals:
  input: expand("15_Trim/LA.intervals.{lg}.trimmed", lg = lg_range)
  output: "16_MareyMapsTrimmed/data.marey.trimmed.gz"
  message: "Concatenating trimmed intervals to {output}"
  shell: "cat {input} | gzip -c > {output}"


rule plot_trimmedintervals:
  input: "16_MareyMapsTrimmed/data.marey.trimmed.gz"
  output: report("16_MareyMapsTrimmed/LepAnchor.mareymaps.pdf", category = "Trimmed Marey Maps")
  message: "Plotting results of edge trimming"
  shell: "LASummary.r {input}"