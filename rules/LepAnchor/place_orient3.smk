
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
    gunzip -fc {input.chain} | java -cp $CONDA_PREFIX/bin/lepanchor PlaceAndOrientContigs chromosome={params.chrom} numThreads={threads} $(awk -f $CONDA_PREFIX/bin/pickorientation.awk {input.chrom}) bed={input.bedfile} map={input.lift} chain=- paf={input.paf} proximity={input.prox} evaluateAnchoring={input.propogated} improveAnchoring=1 {params.datatype} {params.extras} > {output.chrom} 2> {output.errors}
    """
    
rule prune_contigblocks:
  input: "10_PlaceAndOrientContigs/3_orient/chr.{lg_range}.la"
  output: 
    chrom = "10_PlaceAndOrientContigs/pruned/chr.{lg_range}.pruned.la",
    err = "10_PlaceAndOrientContigs/pruned/err/chr.{lg_range}.pruned.err"
  message: "Pruning contig blocks without map support and removing overlaps"
  params:
    chrom = lg
  shell: "awk -f $CONDA_PREFIX/bin/prune.awk {input} > {output.chrom} 2> {output.err}"

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
    awk -f $CONDA_PREFIX/bin/removeOverlaps.awk {input.bedfile} {input.prunedchrom} > {output.overlaps}
    """