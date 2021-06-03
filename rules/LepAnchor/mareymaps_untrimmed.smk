rule mareymap_data:
  input:
    lift = "10_PlaceAndOrientContigs/liftover.la",
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range)
  output: 
    mareydata = "13_MareyMapsUntrimmed/data.marey.gz",
    sexavg = "13_MareyMapsUntrimmed/data.marey.sexavg.gz"
  log: report("13_MareyMapsUntrimmed/missing_scaffolds.txt", category = "Logs")
  message: 
    """
    Marey map interval data
    ===========================================================
    first points in uncertainty intervals  | {output.mareydata}
    midpoints in uncertainty intervals     | {output.sexavg}  
    """
  params:
    chrom = lg
  shell:
    """
    for c in $(seq 1 {params.chrom})
    do
      awk -vn=$c '($3==n)' {input.lift} | awk -f software/LepAnchor/scripts/liftover.awk 11_AGP/contigs/chr.$c.agp - | awk -vm=1 '(/LG/ && NF>=4){{if (NF==4) $5=$4;print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$5}}' | gzip
    done > {output.mareydata} 2> {log}
   
    for c in $(seq 1 {params.chrom})
    do
      awk -vn=$c '($3==n)' {input.lift} | awk -f software/LepAnchor/scripts/liftover.awk 11_AGP/contigs/chr.$c.agp - | awk -vm=1 '(/LG/ && NR>=4){{if (NF>4) s=0.5; else s=1;print $1"\t"$2"\t"$3"\t"m"\t"s*($4+$5)}}' | gzip
    done > {output.sexavg} 2> /dev/null
    """


rule mareymaps:
  input:
    data = "13_MareyMapsUntrimmed/data.marey.gz",
    sexavg = "13_MareyMapsUntrimmed/data.marey.sexavg.gz",
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range)
  output: 
    indiv_plots = report(expand("13_MareyMapsUntrimmed/LG.{lgs}.mareymap.png", lgs = lg_range), category = "Marey Maps"),
    summary = report("13_MareyMapsUntrimmed/LepAnchor.mareymaps.pdf", category = "Marey Maps") ,
    sequential = report("13_MareyMapsUntrimmed/LepAnchor.sequentialmaps.pdf", category = "Sequential Maps"),
    SAsummary = report("13_MareyMapsUntrimmed/LepAnchor.sexavg.mareymaps.pdf", category = "Marey Maps Sex Avg"),
    SAsequential = report("13_MareyMapsUntrimmed/LepAnchor.sexavg.sequentialmaps.pdf", category = "Sequential Maps Sex Avg")
  message: "Creating Marey Maps"
  shell: 
    """
    Rscript software/LepAnchor/scripts/plot_marey.R {input.data} 11_AGP/contigs
    Rscript scripts/LASummary.r {input.data} true
    Rscript scripts/LASummarySexAvg.r {input.sexavg}
    """