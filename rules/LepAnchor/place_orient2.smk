rule place_orient2:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_PlaceAndOrientContigs/map.propogated.bed",
    paf = paf,
    prox = proximity,
    lift = "10_PlaceAndOrientContigs/liftover.la",
    chrom = "10_PlaceAndOrientContigs/1_orient/chr.{lg_range}.la"
  output:
    chrom = "10_PlaceAndOrientContigs/2_orient/chr.{lg_range}.la",
    chromerr = "10_PlaceAndOrientContigs/2_orient/logs/chr.{lg_range}.err"
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra,
    datatype = data_type
  threads: 2
  message: "Running 2nd round of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp $CONDA_PREFIX/bin/lepanchor PlaceAndOrientContigs chromosome={params.chrom} numThreads={threads} $(awk -f $CONDA_PREFIX/bin/pickorientation.awk {input.chrom}) bed={input.bedfile} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} > {output.chrom} 2> {output.chromerr}
    """

rule propogate2:
  input:
    placed = expand("10_PlaceAndOrientContigs/2_orient/chr.{lgs}.la", lgs = lg_range),
    bedfile = "10_PlaceAndOrientContigs/map.propogated.bed"
  output:
    prop = expand("10_PlaceAndOrientContigs/propogate/propogated.{lgs}.la", lgs = lg_range),
    propogated = "10_PlaceAndOrientContigs/map.propogated2.bed",
    iter1 = temp("10_PlaceAndOrientContigs/tmp1.la"),
    iter2 = temp("10_PlaceAndOrientContigs/tmp2.la"),
    iter3 = temp("10_PlaceAndOrientContigs/tmp3.la")
  message: "Second round of propogation"
  shell:
    """
    awk -f $CONDA_PREFIX/bin/propagate.awk {input.placed} > {output.iter1}
    awk -f $CONDA_PREFIX/bin/propagate.awk {output.iter1} > {output.iter2}
    i=2

    while ! cmp -s "10_PlaceAndOrientContigs/tmp$i.la" "10_PlaceAndOrientContigs/tmp$(( $i-1 )).la" ;do
      awk -f $CONDA_PREFIX/bin/propagate.awk 10_PlaceAndOrientContigs/tmp$i.la > 10_PlaceAndOrientContigs/tmp$[$i+1].la
      i=$[$i+1]
    done

    #create prop*.la
    awk -f $CONDA_PREFIX/bin/propagate2.awk 10_PlaceAndOrientContigs/tmp$i.la | awk '(/^[^#]/ && NF>=8){{++d[$1"\t"($7+0)"\t"($8+0)]; data[++line]=$0}}END{{for (i=1; i<=line; ++i) {{$0=data[i];if (d[$1"\t"($7+0)"\t"($8+0)] == 1) {{fn="10_PlaceAndOrientContigs/propogate/propogated."$5".la";print $0>fn}}}}}}'
    awk '{{print $1"\t"($7+0)"\t"($8+0)"\t?\t"$5}}' {output.prop} | awk -f $CONDA_PREFIX/bin/pickbed.awk - {input.bedfile} > {output.propogated}
    """