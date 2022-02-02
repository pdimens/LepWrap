#rule liftover2:
#  input:
#    chain = "9_Chain/chainfile.gz",
#    intervals = "10_PlaceAndOrientContigs/lepanchor.input",
#    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes.before",
#    haplos2 = "10_PlaceAndOrientContigs/suspected.haplotypes.after",
#    lengths = "10_PlaceAndOrientContigs/contigs.length"
#  output:
#    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes.all",
#    lift = report("10_PlaceAndOrientContigs/liftover.nohaplotypes.la", category = "Lifted Intervals"),
#    sortedlift = report("10_PlaceAndOrientContigs/liftover.sorted.nohaplotypes.la", category = "Lifted Intervals"),
#    mapfile = "10_PlaceAndOrientContigs/map.nohaplotypes.clean",
#    bedfile = "10_PlaceAndOrientContigs/map.nohaplotypes.bed",
#    unused = "10_PlaceAndOrientContigs/not_used.nohaplotypes.txt",
#    chr0 = "10_PlaceAndOrientContigs/chr0.nohaplotypes.bed",
#    mapextra = "10_PlaceAndOrientContigs/map.nohaplotypes.extra.bed"
#  message: "Recreating bedfile omitting haplotypes discovered from PlaceAndOrientContigs"
#  params:
#    chrom = lg
#  shell:
#    """
#    cat {input.haplos} {input.haplos2} > {output.haplos}
#    gunzip -fc {input.chain} | java -cp software/LepAnchor LiftoverHaplotypes map={input.intervals} haplotypes={output.haplos} chain=- > {output.lift}
#    cat {output.lift} | sort -V -k 1,1 -k 2,2n > {output.sortedlift}
#    java -cp software/LepAnchor CleanMap map={output.sortedlift} > {output.mapfile}
#    java -cp software/LepAnchor Map2Bed map={output.mapfile} contigLength={input.lengths} > {output.bedfile}
#    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {output.haplos}; cut -f 1 {output.bedfile}) > {output.unused}
#    grep -w -F -f {output.unused} {input.lengths} | awk -vn={params.chrom} '{{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=n;++i) print s i}}' > {output.chr0}
#    cat {output.bedfile} {output.chr0} > {output.mapextra}
#    """


rule place_orient2:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map.propogated.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.la",
    chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la"
  output:
    chrom = "10_PlaceAndOrientContigs/orient_2/chr.{lg_range}.la",
    haplos = "10_PlaceAndOrientContigs/orient_2/haplotypes/chr.{lg_range}.haplo.suspected",
    chromerr = report("10_PlaceAndOrientContigs/orient_2/logs/chr.{lg_range}.la.log", category = "Anchoring II Logs")
  log:
    haplos = "10_PlaceAndOrientContigs/orient_2/haplotypes/chr.{lg_range}.haplo.all",
    errors = "10_PlaceAndOrientContigs/orient_2/errors/chr.{lg_range}.errors"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type,
    haplo = haplo_limit
  threads: 3
  message: "Running 2nd round of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs numThreads={threads} $(awk -f software/LepAnchor/scripts/pickorientation.awk {input.chrom}) bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} > {output.chrom} 2> {output.chromerr}
    sort -n -r {output.chromrr} | awk '($NF=="haplotype" && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}' > {log.haplos}
    sort -n -r {output.chromerr} | awk -vlimit={params.haplo} '($NF=="haplotype" && ($1>=($4-$3+1-limit)/limit) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){{h[$2,$3,$4]; print}}'  > {output.haplos}  
    grep "error$" {log.chromerr} > {log.errors}
    """

  rule propogate2:
    input:
      placed = expand("10_PlaceAndOrientContigs/orient_2/logs/chr.{lgs}.la.log", lgs = lg_range),
      bedfile = "10_PlaceAndOrientContigs/map.propogated.bed"
    output:
      prop = expand("10_PlaceAndOrientContigs/propogate/propogated.{lgs}.la", lgs = range(lg + 1)),
      propogated = "10_PlaceAndOrientContigs/map.propogated2.bed",
      iter1 = temp("10_PlaceAndOrientContigs/tmp1.la"),
      iter2 = temp("10_PlaceAndOrientContigs/tmp2.la")
    message: "Second round of propogation"
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
      awk -f software/LepAnchor/scripts/propagate2.awk tmp$i.la | awk '(/^[^#]/ && NF>=8){{++d[$1"\t"($7+0)"\t"($8+0)]; data[++line]=$0}}END{{for (i=1; i<=line; ++i) {{$0=data[i];if (d[$1"\t"($7+0)"\t"($8+0)] == 1) {{fn="10_PlaceAndOrientContigs/propogate/propogated."$5".la";print $0>fn}}}}}}'
      awk '{{print $1"\t"($7+0)"\t"($8+0)"\t?\t"$5}}' {output.prop} | awk -f software/LepAnchor/scripts/pickbed.awk - {input.bedfile} > {output.propogated}
      """