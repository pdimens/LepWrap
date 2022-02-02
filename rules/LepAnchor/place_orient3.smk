rule place_orient3:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map.propogated2.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.nohaplotypes.la",
    chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la",
    propogated = "10_PlaceAndOrientContigs/propogate/propogated.{lg_range}.la"
  output:
    chrom = "10_PlaceAndOrientContigs/orient_3/chr.{lg_range}.la",
    haplos = "10_PlaceAndOrientContigs/orient_3/haplotypes/chr.{lg_range}.haplo.suspected"
  log:
    chrom = report("10_PlaceAndOrientContigs/orient_3/logs/chr.{lg_range}.la.log", category = "Anchoring III Logs"),
    haplos = "10_PlaceAndOrientContigs/orient_3/haplotypes/chr.{lg_range}.haplo.all",
    errors = "10_PlaceAndOrientContigs/orient_3/errors/chr.{lg_range}.errors"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type,
    haplo = haplo_limit
  message: "Running 3rd round of PlaceAndOrientContigs for linkage group {params.chrom}"
  threads: 3
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs numThreads={threads} $(awk -f software/LepAnchor/scripts/pickorientation.awk {input.chrom}) bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} evaluateAnchoring={input.propogated} improveAnchoring=1 {params.datatype} {params.extras} > {output.chrom} 2> {log.chrom}
    sort -n -r {log.chrom} | awk '($NF=="haplotype" && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}' > {log.haplos}
    sort -n -r {log.chrom} | awk -vlimit={params.haplo} '($NF=="haplotype" && ($1>=($4-$3+1-limit)/limit) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}'  > {output.haplos}  
    grep "error$" {log.chrom} > {log.errors}
    """


rule prune:
  input: 
    oriented = expand("10_PlaceAndOrientContigs/orient_3/chr.{lgs}.la", lgs = lg_range),
    bedfile = "10_PlaceAndOrientContigs/map.propogated2.bed"
  output: 
    pruned = report("10_PlaceAndOrientContigs/orient_3/pruned.la", category = "Logs"),
    cleaned = report("10_PlaceAndOrientContigs/overlaps.removed.la", category = "Logs")
  message: "Pruning contig blocks without map support and removing overlaps"
  params:
    chrom = lg
  shell:
    """
    for i in $(seq {params.chrom})
    do
      awk -f software/LepAnchor/scripts/prune.awk 10_PlaceAndOrientContigs/orient_3/ichr.$i.la > 10_PlaceAndOrientContigs/orient_3/ichr.${{i}}.pruned.la
    done 2> {output.pruned}
    awk -f software/LepAnchor/scripts/removeOverlaps.awk {input.bedfile} 10_PlaceAndOrientContigs/orient_3/ichr.*.pruned.la > {output.cleaned}
    """