rule place_orient:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map_extra.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.la"
  output:
    chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la",
    haplos = "10_PlaceAndOrientContigs/orient_1/haplotypes/chr.{lg_range}.haplo.suspected",
    chromerr = report("10_PlaceAndOrientContigs/orient_1/logs/chr.{lg_range}.la.log", category = "Anchoring I Logs")
  log:
    errors = "10_PlaceAndOrientContigs/orient_1/errors/chr.{lg_range}.errors",
    haplos = "10_PlaceAndOrientContigs/orient_1/haplotypes/chr.{lg_range}.haplo.all"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type,
    haplo = haplo_limit
  threads: 3
  message: "Running the 1st round of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs numThreads={threads} bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} > {output.chrom} 2> {output.chromerr}
    sort -n -r {log.chrom} | awk '($NF=="haplotype" && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}' > {log.haplos}
    sort -n -r {log.chrom} | awk -vlimit={params.haplo} '($NF=="haplotype" && ($1>=($4-$3+1-limit)/limit) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}'  > {output.haplos}  
    grep "error$" {log.chrom} > {log.errors}
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
    #awk -f propagate4.awk pass=1 chr*[0-9].la pass=2 chr*[0-9].la.err|awk -f pickbed.awk - map_extra.bed >map_extra_prop.bed
    """

rule mergehaplos:
  input: expand("10_PlaceAndOrientContigs/orient_1/haplotypes/chr.{lg}.haplo.suspected", lg = lg_range)
  output: "10_PlaceAndOrientContigs/suspected.haplotypes.after"
  message: "Merging suspected haplotype contig information from the linkage groups"
  shell: "cat {input} | sort | uniq > {output}"
