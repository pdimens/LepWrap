rule place_orient:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map_extra.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.la"
  output:
    chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la",
    chromerr = report("10_PlaceAndOrientContigs/orient_1/logs/chr.{lg_range}.la.log", category = "Anchoring I Logs")
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type
  threads: 2
  message: "Running the 1st round of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs chromosome={params.chrom} numThreads={threads} bed={input.bedfile} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} > {output.chrom} 2> {output.chromerr}
    """

rule propogate1:
  input: 
    placed = expand("10_PlaceAndOrientContigs/orient_1/chr.{lgs}.la", lgs = lg_range),
    errors = expand("10_PlaceAndOrientContigs/orient_1/logs/chr.{lgs}.la.log", lgs = lg_range),
    bedfile = "10_PlaceAndOrientContigs/map_extra.bed"
  output:
    propogated = "10_PlaceAndOrientContigs/map.propogated.bed",
  message: "First round of propogation with propogate4.awk"
  shell:
    """
    awk -f software/LepAnchor/scripts/propagate4.awk pass=1 {input.placed} pass=2 {input.errors} | awk -f software/LepAnchor/scripts/pickbed.awk - {input.bedfile} > {output.propogated}
    """
