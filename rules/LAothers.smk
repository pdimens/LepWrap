if [[ ! $REF =~ ^$ ]];then
	echo "calculating contigs.length file..."
	echo "gunzip -fc $REF | awk -f LA/scripts/contigLength.awk > contigs.length" | bash
fi



rule contiglengths:
  input: geno
  output: "10_Anchoring/contigs.length"
  message: "Getting contig lengths"
  shell: "gunzip -fc {input} | awk -f LA/scripts/contigLength.awk > {output}"


echo "finding full haplotypes..."
echo "gunzip -fc $CHAIN | awk -f scripts/findFullHaplotypes.awk > fullHaplotypes50.txt" | bash
wc -l fullHaplotypes50.txt

rule find_haplotypes:
  input: "9_Chain/chainfile.gz"
  output: "10_Anchoring/fullHaplotypes50.txt"
  message: "Finding full haplotypes (potential chimeric contigs)"
  shell: 
    """
    gunzip -fc {input} | awk -f LA/scripts/findFullHaplotypes.awk > {output}
    echo "Detected $(wc -l {output}) potential chimeric contigs
    """

rule find_haplotypes:
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
    grep -w -F -f 10_Anchoring/not_used.txt {input.lengths} | awk -vn=$CHR '{{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=n;++i) print s i}}' > chr0.bed
    cat {input.bedfile} chr0.bed > {output.bedfile}
    """

rule place_orient:
  input:
    chain = "9_Chain/chainfile.gz",
    bedfile = "10_Anchoring/map_extra.bed",
    paf = paf,
    prox = proximity,
    lift = "10_Anchoring/liftover.la"
  output:
    chrom = "10_Anchoring/chr.{lg_range}.la"
  log:
    chrom = "10_Anchoring/logs/chr.{lg_range}.la.err"
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
    awk '(NR==FNR){{print;c[$1]}}(NR!=FNR && !($1 in c)){{print $1 "\t" $7+0 "\t" $8+0"\t?\t"$5}}' {input.bedfile} propogated.[1-9]*.la > {output.propogated}
    rm 10_Anchoring/tmp*.la
    """


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