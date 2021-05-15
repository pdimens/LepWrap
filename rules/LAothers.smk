rule contiglengths:
  input: geno
  output: "10_Anchoring/contigs.length"
  message: "Getting contig lengths"
  shell: "gunzip -fc {input} | awk -f LA/scripts/contigLength.awk > {output}"

rule find_haplotypes:
  input: "9_Chain/chainfile.gz"
  output: "10_Anchoring/fullHaplotypes50.txt"
  message: "Finding full haplotypes (potential chimeric contigs)"
  shell: 
    """
    gunzip -fc {input} | awk -f LA/scripts/findFullHaplotypes.awk > {output}
    echo "Detected $(wc -l {output}) potential chimeric contigs
    """

rule liftover:
  input: 
    chain = "9_Chain/chainfile.gz",
    intervals = "10_Anchoring/lepmap3_intervals.la",
    haplos = "10_Anchoring/fullHaplotypes50.txt"
  output: 
    lift = "10_Anchoring/liftover.la",
    sortedlift = "10_Anchoring/liftover.sorted.la"
  message: "Running liftoverHaplotypes for the input maps"
  shell: 
    """
    gunzip -fc {input.chain} | java -cp LA LiftoverHaplotypes map={input.intervals} haplotypes={input.haplos} chain=- > {output.lift}
    cat {output.lift} | sort -V -k 1,1 -k 2,2n > {sortedlift}
    """

rule cleanmap:
  input: "10_Anchoring/liftover.sorted.la"
  output: "10_Anchoring/map_all.clean"
  message: "Running CleanMap"
  shell: "java -cp LA CleanMap map={input} > {output}"

rule map2bed:
  input: 
    cleanmap = "10_Anchoring/map_all.clean",
    lengths = "10_Anchoring/contigs.length",
  output: "10_Anchoring/map.bed"
  message: "Running Map2Bed"
  shell: "java -cp LA Map2Bed map={input.cleanmap} contigLength={input.lengths} > {output}"

rule ungrouped:
  input:
    lengths = "10_Anchoring/contigs.length",
    haplos = "10_Anchoring/fullHaplotypes50.txt",
    bedfile = "10_Anchoring/map.bed"
  output:
    bedfile = "10_Anchoring/map_extra.bed"
  message: "Finding contigs not put into chromosomes"
  shell:
    """
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {input.haplos}; cut -f 1 {input.bedfile}) > 10_Anchoring/not_used.txt
    grep -w -F -f 10_Anchoring/not_used.txt {input.lengths} | awk -vn=$CHR '{{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=n;++i) print s i}}' > 10_Anchoring/chr0.bed
    cat {input.bedfile} 10_Anchoring/chr0.bed > {output.bedfile}
    """

rule place_orient:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_Anchoring/map_extra.bed",
    paf = paf,
    prox = proximity,
    lift = "10_Anchoring/liftover.la"
  output:
    chrom = "10_Anchoring/orient_1/chr.{lg_range}.la"
  log:
    chrom = "10_Anchoring/orient_1/logs/chr.{lg_range}.la.err"
  params:
    chrom = "{lg_range}"
  message: "Running PlaceAndOrientContigs"
  shell:
    """
    echo "gunzip -fc {input.chain} | java -cp LA PlaceAndOrientContigs bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} keepEmptyIntervals=1 > {output} 2> {log}"
    """

rule propogate:
  input: 
    placed = expand("10_Anchoring/chr.{lgs}.la", lgs = lg_range),
    bedfile = "10_Anchoring/map.bed"
  output:
    propogated = "10_Anchoring/map_propogated.bed",
    tmp_prop = temp(expand("10_Anchoring/propogated.{lgs}.la", lgs = lg_range))
  message: "Propogating ...something"
  shell:
    """
    awk -f LA/scripts/propagate.awk chr*.la > 10_Anchoring/tmp1.la
    awk -f LA/scripts/propagate.awk tmp1.la > 10_Anchoring/tmp2.la
    i=2

    while ! cmp -s "10_Anchoring/tmp$i.la" "10_Anchoring/tmp$(( $i-1 )).la" ;do
	    awk -f LA/scripts/propagate.awk 10_Anchoring/tmp$i.la > 10_Anchoring/tmp$[$i+1].la
	    i=$[$i+1]
    done
    #create prop*.la
    awk '/^[^#]/{{++d[$1 "\t" $7+0 "\t" $8+0]; data[++line]=$0}}END{{for (i = 1; i <= line; ++i) {{$0=data[i];if (d[$1 "\t" $7+0 "\t" $8+0] == 1) fn="10_Anchoring/propogated."$5".la"; else if ($5==1) fn="10_Anchoring/propogated.0.la"; else fn=""; if (fn != "") print $0>fn}}' 10_Anchoring/tmp$i.la

    #create a new bed by combining propogated.[1-9]*.la and map.bed
    awk '(NR==FNR){{print;c[$1]}}(NR!=FNR && !($1 in c)){{print $1 "\t" $7+0 "\t" $8+0"\t?\t"$5}}' {input.bedfile} {output.tmp_prop} > {output.propogated}
    rm 10_Anchoring/tmp*.la
    """

rule place_orient2:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_Anchoring/map_propogated.bed",
    paf = paf,
    prox = proximity,
    lift = "10_Anchoring/liftover.la"
  output:
    chrom = "10_Anchoring/orient_2/ichr.{lg_range}.la"
  log:
    chrom = "10_Anchoring/orient_2/logs/ichr.{lg_range}.la.err"
  params:
    chrom = "{lg_range}"
  message: "Running a second iteration of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    echo "gunzip -fc {input.chain} | java -cp LA PlaceAndOrientContigs bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} keepEmptyIntervals=1 > {output.chrom} 2> {log.chrom}
    """

rule prune:
  input: 
    oriented = expand("10_Anchoring/orient_2/ichr.{lgs}.la", lgs = lg_range),
    bedfile = "10_Anchoring/map_propogated.bed"
  output: 
    pruned = "10_Anchoring/orient_2/pruned.la",
    cleaned = "10_Anchoring/overlaps_rm.la"
  message: "Pruning contig blocks without map support and removing overlaps"
  params:
    chrom = lg
  shell:
    """
    for i in $(seq {params.chrom})
    do
      awk -f LA/scripts/prune.awk 10_Anchoring/orient_2/ichr.$i.la > 10_Anchoring/orient_2/ichr.${i}.pruned.la
    done 2> {output.pruned}
    awk -f LA/scripts/removeOverlaps.awk {input.bedfile} 10_Anchoring/orient_2/ichr.${i}.pruned.la > {output.cleaned}
    """

rule construct_agp:
  input:
    cleaned = "10_Anchoring/overlaps_rm.la"
  output:
    agp = "10_Anchoring/agp/chr.{lg_range}.agp",
    scaff_agp = "10_Anchoring/agp_scaffolds/chr.{lg_range}.scaffolds.agp"
  message: "Creating AGP files for linkage group {lg_range}"
  params:
    chrom = "{lg_range}"
  shell:
    """
    awk -vn={params.chrom} '($5==n)' {input.cleaned} | awk -vprefix="LG" -vlg={params.chrom} -f LA/scripts/makeagp_full2.awk - > 10_Anchoring/agp/chr.{params.chrom}.agp
    awk -vn={params.chrom} '($5==n)' {input.cleaned} | awk -vprefix="LG" -vlg={params.chrom} -f LA/scripts/makeagp2.awk - > 10_Anchoring/agp_scaffolds/chr.{params.chrom}.scaffolds.agp
    """

rule unused:
  input:
    lengths = "10_Anchoring/contigs.length",
    haplos = "10_Anchoring/fullHaplotypes50.txt",
    agp = expand("10_Anchoring/agp/chr.{lgs}.agp", lgs = lg_range),
    scaff_agp = expand("10_Anchoring/agp_scaffolds/chr.{lgs}.scaffolds.agp", lgs = lg_range)
  output: 
    txt = "10_Anchoring/not_used_final.txt",
    agp = "10_Anchoring/not_used.agp",
    final_agp = "10_Anchoring/REF_LA.agp",
    scaff_agp = "10_Anchoring/REF_LA_scaffolds.agp"
  message: "Finding unused contigs"
  shell:
    """
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {inputs.haplos};awk '($5!="U"){print $6}' {input.agp}) > {output.txt}
    grep -F -w -f {output.txt} {inputs.lengths} | awk '{print $1,1,$2,1,"W",$1,1,$2,"+"}' > {output.agp}
    cat {input.agp} {output.agp} > {output.final_agp}
    cat {input.scaff_agp} {output.agp} > {output.scaff_agp}
    """


rule build_fasta:
  input:
    assembly = geno
    agp = "10_Anchoring/REF_LA.agp",
    scaff_agp = "10_Anchoring/REF_LA_scaffolds.agp"
  output:
    fasta = "10_Anchoring/Anchored.contigs.fa.gz",
    scaff = "10_Anchoring/Anchored.scaffolds.fa.gz"
  message: "Constructing final fasta files"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f LA/scripts/makefasta.awk - {input.agp} | gzip > {output.fasta}
    gunzip -fc {input.assembly} | awk -f LA/scripts/makefasta.awk - {input.scaff_agp} | gzip > {output.scaff}
    """

rule mareymaps:
  input:
  output: "11_MareyMaps/"
  message: "Creating Marey map plots"
  params:
    chrom = lg
  shell:
    """
    j=1
    for m in $MAP
    for c in $(seq {params.chrom})
    do
    awk -vn=$c '($3==n)' $m.liftover | awk -f LA/scripts/liftover.awk 10_Anchoring/chr$c.agp - |awk -vm=$j '(/LG/ && NF>=4){if (NF==4) $5=$4;print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$5}' | gzip > marey.data.gz
    done
    Rscript LA/scripts/plot_marey.R
    """

    j=1
    for m in $MAP
    do
      for c in $(seq $CHR)
      do
      awk -vn=$c '($3==n)' $m.liftover | awk -f LA/scripts/liftover.awk chr$c.agp - |awk -vm=$j '(/LG/ && NF>=4){if (NF==4) $5=$4;print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$5}'
      done
      j=$[$j + 1]
    done | gzip > marey.data.gz

    Rscript scripts/plot_marey.R

rule lepanchor:
  input:
    genome = {assembly},
    paf_file = {paf},
    chain = "",
    intervals = "10_Anchoring/lepmap3_intervals.la"
  output:
  log: "10_Anchoring/LepAnchor.log"
  message: "Running LepAnchor"
  threads: 30
  shell: 
    """
    if [ {input.paf_file} == "none" ]; then
      lepanchor_wrapper.sh -t {threads} -f {input.genome} -n {lg} -c {input.chain} -m {input.intervals} 2> {log}
    else
      lepanchor_wrapper.sh -t {threads} -f {input.genome} -n {lg} -c {input.chain} -m {input.intervals} -p {input.paf} 2> {log}
    fi
    """