rule place_orient3:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map.propogated2.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.la",
    chrom = "10_PlaceAndOrientContigs/1_orient/chr.{lg_range}.la",
    propogated = "10_PlaceAndOrientContigs/propogate/propogated.{lg_range}.la"
  output:
    chrom = "10_PlaceAndOrientContigs/3_orient/chr.{lg_range}.la",
    errors = "10_PlaceAndOrientContigs/3_orient/errors/chr.{lg_range}.err"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type
  message: "Running 3rd round of PlaceAndOrientContigs for linkage group {params.chrom}"
  threads: 2
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs chromosome={params.chrom} numThreads={threads} $(awk -f software/LepAnchor/scripts/pickorientation.awk {input.chrom}) bed={input.bedfile} map={input.lift} chain=- paf={input.paf} proximity={input.prox} evaluateAnchoring={input.propogated} improveAnchoring=1 {params.datatype} {params.extras} > {output.chrom} 2> {output.errors}
    """
    
rule prune_contigblocks:
  input: "expand(10_PlaceAndOrientContigs/3_orient/chr.{lg_range}.la")
  output: 
    chrom = "10_PlaceAndOrientContigs/pruned/chr.{lg_range}.pruned.la",
    err = "10_PlaceAndOrientContigs/pruned/err/chr.{lg_range}.pruned.err"
  message: "Pruning contig blocks without map support and removing overlaps"
  params:
    chrom = lg
  shell: "awk -f software/LepAnchor/scripts/prune.awk {input} > {output.chrom} 2> {output.err}"

rule prune_post:
  input:
    bedfile = "10_PlaceAndOrientContigs/map.propogated2.bed",
    prunedchrom = expand("10_PlaceAndOrientContigs/pruned/chr.{lgs}.pruned.la", lgs = lg_range),
    prunederr = expand("10_PlaceAndOrientContigs/pruned/err/chr.{lgs}.pruned.err", lgs = lg_range)
  output: 
    overlaps = "10_PlaceAndOrientContigs/overlaps.removed.la",
    pruned = "10_PlaceAndOrientContigs/pruned.la"
  message: "Removing overlaps"
  threads: 1
  shell:
    """
    cat {input.prunederr} > {output.pruned}
    awk -f software/LepAnchor/scripts/removeOverlaps.awk {input.bedfile} {input.prunedchrom} > {output.overlaps}
      """

#rule find_haplotypes2:
#  input:
#    errors = expand("10_PlaceAndOrientContigs/3_orient/errors/chr.{lgs}.err", lgs = lg_range)
#  output:
#    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes"
#  message: "Finding full haplotypes"
#  params:
#    haplo = haplo_limit
#  shell:
#    """
#    awk '($NF=="haplotype")' {input.errors} | 
#      sort -n -r | 
#      awk -vlimit={params.haplo} '($NF=="haplotype" && ($1>=($4-$3+1-limit)/limit) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}' > {output.haplos}
#    """
#
#rule liftoverHaplotypes:
#  input:
#    chain = "9_Chain/chainfile.gz",
#    chrom = "10_PlaceAndOrientContigs/1_orient/chr.{lg_range}.la",
#    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes"
#  output: "10_PlaceAndOrientContigs/liftover/chr.{lg_range}.liftover"
#  message: "Running liftoverHaplotypes for {input.chrom}"
#  threads: 1
#  shell:
#    """
#    gunzip -fc {input.chain} | java -cp software/LepAnchor LiftoverHaplotypes map={input.chrom} haplotypes={input.haplos} chain=- > {output}
#    """
#
#rule removehaplotypes:
#  input:
#    mapfile = "10_PlaceAndOrientContigs/map.propogated2.bed", 
#    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes"
#  output:
#    bedfile = "10_PlaceAndOrientContigs/map.propogated2.nohaplo.bed"
#  message: "Removing haplotypes from the map"
#  threads: 1
#  shell: "awk -f software/LepAnchor/scripts/removeHaplotypes.awk {input} > {output}"
#