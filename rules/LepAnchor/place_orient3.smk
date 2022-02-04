rule place_orient3:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map.propogated2.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.la",
    chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la",
    propogated = "10_PlaceAndOrientContigs/propogate/propogated.{lg_range}.la"
  output:
    chrom = "10_PlaceAndOrientContigs/orient_3/chr.{lg_range}.la",
    errors = "10_PlaceAndOrientContigs/orient_3/errors/chr.{lg_range}.err"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type,
  message: "Running 3rd round of PlaceAndOrientContigs for linkage group {params.chrom}"
  threads: 2
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs chromosome={params.chrom} numThreads={threads} $(awk -f software/LepAnchor/scripts/pickorientation.awk {input.chrom}) bed={input.bedfile} map={input.lift} chain=- paf={input.paf} proximity={input.prox} evaluateAnchoring={input.propogated} improveAnchoring=1 {params.datatype} {params.extras} > {output.chrom} 2> {output.errors}
    """
    
rule find_haplotypes2:
  input:
    errors = expand("10_PlaceAndOrientContigs/orient_3/errors/chr.{lgs}.err", lgs = lg_range)
  output:
    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes"
  message: "Finding full haplotypes"
  params:
    haplo = haplo_limit
  shell:
    """
    awk '($NF=="haplotype")' {input.errors} | 
      sort -n -r | 
      awk -vlimit={params.haplo} '($NF=="haplotype" && ($1>=($4-$3+1-limit)/limit) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}' > {output.haplos}
    """

rule liftoverHaplotypes:
  input:
    chain = "9_Chain/chainfile.gz",
    chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la",
    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes"
  output: "10_PlaceAndOrientContigs/liftover/chr.{lg_range}.liftover"
  message: "Running liftoverHaplotypes for {input.chrom}"
  threads: 1
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor LiftoverHaplotypes map={input.chrom} haplotypes={input.haplos} chain=- > $i.liftover
    """

rule removehaplotypes:
  input:
    mapfile = "10_PlaceAndOrientContigs/map.propogated2.bed", 
    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes"
  output:
    bedfile = "10_PlaceAndOrientContigs/map.propogated2.nohaplo.bed"
  message: "Removing haplotypes from the map"
  threads: 1
  shell: "awk -f software/LepAnchor/scripts/removeHaplotypes.awk {input} > {output}"
