rule construct_agp:
  input:
    cleaned = "10_PlaceAndOrientContigs/overlaps.removed.la"
  output:
    agp = report("11_AGP/contigs/chr.{lg_range}.agp", category = "Contig AGP Files"),
    scaff_agp = report("11_AGP/scaffolds/chr.{lg_range}.scaffolds.agp", category = "Scaffold AGP Files")
  message: "Creating AGP files for linkage group {params.chrom}"
  params:
    chrom = "{lg_range}"
  shell:
    """
    awk -vn={params.chrom} '($5==n)' {input.cleaned} | awk -vprefix="LG" -vlg={params.chrom} -f software/LepAnchor/scripts/makeagp_full2.awk - > {output.agp}
    awk -vn={params.chrom} '($5==n)' {input.cleaned} | awk -vprefix="LG" -vlg={params.chrom} -f software/LepAnchor/scripts/makeagp2.awk - > {output.scaff_agp}
    """


rule unused:
  input:
    lengths = "10_PlaceAndOrientContigs/contigs.length",
    haplos = "10_PlaceAndOrientContigs/fullHaplotypes50.txt",
    cleaned = "10_PlaceAndOrientContigs/overlaps.removed.la",
  output: 
    txt = "11_AGP/not_used_final.txt",
    agp = "11_AGP/not_used.agp"
  message: "Finding unused contigs"
  shell:
    """
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f2 {input.haplos}; cut -f1 {input.cleaned}) > {output.txt}
    grep -F -w -f {output.txt} {input.lengths} | awk '{{print $1,1,$2,1,"W",$1,1,$2,"+"}}' > {output.agp}
    """


rule build_final_agp:
  input:
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range),
    scaff_agp = expand("11_AGP/scaffolds/chr.{lgs}.scaffolds.agp", lgs = lg_range),
    unused = "11_AGP/not_used.agp",
  output:
    contig_agp = "11_AGP/lepanchor.contigs.only.agp",
    scaff_agp = "11_AGP/lepanchor.scaffolds.only.agp",
    contig_all_agp = "11_AGP/lepanchor.contigs.all.agp",
    scaff_all_agp = "11_AGP/lepanchor.scaffolds.all.agp"
  message: "Generating final AGP files"
  shell:
    """
    cat {input.agp} > {output.contig_agp}
    cat {input.scaff_agp} > {output.scaff_agp}
    cat {input.agp} {input.unused} > {output.contig_all_agp}
    cat {input.scaff_agp} {input.unused} > {output.scaff_all_agp}
    """