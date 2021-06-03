rule propogate:
  input: 
    placed = expand("10_PlaceAndOrientContigs/orient_2/chr.{lgs}.la", lgs = lg_range),
    bedfile = "10_PlaceAndOrientContigs/map.nohaplotypes.bed"
  output:
    propogated = "10_PlaceAndOrientContigs/map.propogated.bed",
    tmp_prop = temp(expand("10_PlaceAndOrientContigs/propogate/propogated.{lgs}.la", lgs = range(lg + 1))),
  message: "Propogating ...something"
  shell:
    """
    awk -f software/LepAnchor/scripts/propagate.awk {input.placed} > 10_PlaceAndOrientContigs/tmp1.la
    awk -f software/LepAnchor/scripts/propagate.awk 10_PlaceAndOrientContigs/tmp1.la > 10_PlaceAndOrientContigs/tmp2.la
    i=2

    while ! cmp -s "10_PlaceAndOrientContigs/tmp$i.la" "10_PlaceAndOrientContigs/tmp$(( $i-1 )).la" ;do
	    awk -f software/LepAnchor/scripts/propagate.awk 10_PlaceAndOrientContigs/tmp$i.la > 10_PlaceAndOrientContigs/tmp$[$i+1].la
	    i=$[$i+1]
    done
    #create prop*.la
    awk '/^[^#]/{{++d[$1 "\t" $7+0 "\t" $8+0]; data[++line]=$0}}END{{for (i = 1; i <= line; ++i) {{$0=data[i];if (d[$1 "\t" $7+0 "\t" $8+0] == 1) fn="10_PlaceAndOrientContigs/propogate/propogated."$5".la"; else if ($5==1) fn="10_PlaceAndOrientContigs/propogate/propogated.0.la"; else fn=""; if (fn != "") print $0>fn}}}}' 10_PlaceAndOrientContigs/tmp$i.la

    #create a new bed by combining propogated.[1-9]*.la and map.nohaplotypes.bed
    awk '(NR==FNR){{print;c[$1]}}(NR!=FNR && !($1 in c)){{print $1 "\t" $7+0 "\t" $8+0"\t?\t"$5}}' {input.bedfile} {output.tmp_prop} > {output.propogated}
    rm 10_PlaceAndOrientContigs/tmp*.la
    """


rule place_orient3:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map.propogated.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.nohaplotypes.la"
  output:
    chrom = "10_PlaceAndOrientContigs/orient_3/ichr.{lg_range}.la",
    haplos = "10_PlaceAndOrientContigs/orient_3/haplotypes/chr.{lg_range}.haplo.suspected"
  log:
    chrom = report("10_PlaceAndOrientContigs/orient_3/logs/ichr.{lg_range}.la.log", category = "Anchoring III Logs"),
    haplos = "10_PlaceAndOrientContigs/orient_3/haplotypes/chr.{lg_range}.haplo.all",
    errors = "10_PlaceAndOrientContigs/orient_3/errors/chr.{lg_range}.errors"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type,
    haplo = haplo_limit
  message: "Running 3rd round of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} > {output.chrom} 2> {log.chrom}
    sort -n -r {log.chrom} | awk '($NF=="haplotype" && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}' > {log.haplos}
    sort -n -r {log.chrom} | awk -vlimit={params.haplo} '($NF=="haplotype" && ($1>=($4-$3+1-limit)/limit) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}'  > {output.haplos}  
    grep "error$" {log.chrom} > {log.errors}
    """


rule prune:
  input: 
    oriented = expand("10_PlaceAndOrientContigs/orient_3/ichr.{lgs}.la", lgs = lg_range),
    bedfile = "10_PlaceAndOrientContigs/map.propogated.bed"
  output: 
    pruned = report("10_PlaceAndOrientContigs/orient_3/pruned.la", category = "Logs"),
    cleaned = report("10_PlaceAndOrientContigs/overlaps_rm.la", category = "Logs")
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