from os import path
import glob

configfile: "config.yaml"

geno = config["assembly"]
paf = config["PAF_file"]
proximity = config["proximity_file"]
lg = config["lg_count"]
lg_range = list(range(1,lg+1))
os_name = config["OS_info"]

rule all:
  input:
    fasta = "10_Anchoring/Anchored.contigs.fa.gz",
    scaff = "10_Anchoring/Anchored.scaffolds.fa.gz",
    mareydata = "11_MareyMaps/marey.data.gz",
    mareymaps = expand("11_MareyMaps/LG.{lgs}.mareymap.png", lgs = lg_range)
  message: "Lep-Anchor has finished. Good luck with the rest of your analyses!"

rule repeatmask:
  input: geno
  output: "8_Repeatmask/repeatmasked.fa.gz"
  log: "8_Repeatmask/Red.log"
  message: "Using Red to repeat-mask the genome assembly {input}"
  threads: 30
  shell:
    """
    file={input}
    if [ "${{file: -3}}" == ".gz" ]; then
      echo "- Assembly is compressed, creating decompressed copy"
      file=$(basename $file .gz)
      gunzip --stdout {input} > $file
    fi
    ext=$(echo $file | rev | cut -d"." -f1 | rev)
    if [ $ext != "fa" ]; then
      echo "- Assembly extension must end in .fa for Red, creating a corrected symlink"
      ln -srf $file ${{file}}.fa
    fi
    echo "- Running Red"
    LA/deps/Red -gnm . -msk 8_Repeatmask -sco 8_Repeatmask -cnd 8_Repeatmask -rpt 8_Repeatmask > {log} 2>> {log}
    echo "- Compressing repeat-masked genome from Red"
    gzip --stdout 8_Repeatmask/*.msk > {output} && rm 8_Repeatmask/*.msk
    """
    
  
rule chain_1:
  input: 
    geno = "8_Repeatmask/repeatmasked.fa.gz",
    ctrl = "LA/deps/all_lastz.ctl",
    scoremtx = "LA/deps/scoreMatrix.q"
  output: 
    out1 = "9_Chain/repeatmaskedx.sizes",
    out2 = "9_Chain/repeatmasked.sizes"
  message: "Running Lastz via HaploMerger2"
  threads: 30
  params:
    os = os_name
  shell:
    """
    OS=$(echo {params} | tr '[:upper:]' '[:lower:]')
    echo "Using the $OS lastz/chainNet binaries"
    if [ $OS == "ubuntu" ]
    then
        export PATH="$PATH:LA/deps/ubuntu"
    elif [ $OS == "centos5" ]
    then
        export PATH="$PATH:LA/deps/centOS5"
    elif [ $OS == "centos6" ]
    then
        export PATH="$PATH:LA/deps/centOS6"
    else
        echo "$OS is not recognized as one of Ubuntu, CentOS5, or CentOS6, defaulting to Ubuntu"
        export PATH="$PATH:LA/deps/ubuntu"
    fi
    ln -srf {input} 9_Chain/
    cd 9_Chain
    ../LA/deps/step1.HM2 repeatmasked {threads}
    """

rule chain_2:
  input: 
    f1 = "9_Chain/repeatmaskedx.sizes",
    f2 = "9_Chain/repeatmasked.sizes"
  output: 
    original = "9_Chain/repeatmasked.repeatmaskedx.result/all.chain.gz",
    slink = "9_Chain/chainfile.gz"
  message: "Running HaploMerger2 to generate the chain file"
  threads: 30
  params:
    os = os_name
  shell:
    """
    OS=$(echo {params} | tr '[:upper:]' '[:lower:]')    echo "Using the $OS lastz/ chainNet binaries"
    echo "Using the $OS lastz/chainNet binaries"
    if [ $OS == "ubuntu" ]
    then
        export PATH="$PATH:LA/deps/ubuntu"
    elif [ $OS == "centos5" ]
    then
        export PATH="$PATH:LA/deps/centOS5"
    elif [ $OS == "centos6" ]
    then
        export PATH="$PATH:LA/deps/centOS6"
    else
        echo "$OS is not recognized as one of Ubuntu, CentOS5, or CentOS6, defaulting to Ubuntu"
        export PATH="$PATH:LA/deps/ubuntu"
    fi
    
    cd 9_Chain
    ../LA/deps/step2.HM2 repeatmasked {threads} && rm -r repeatmasked.repeatmaskedx.result/raw.axt
    ln -sr {output.original} {output.slink}
    """

rule extract_markers:
    input: "2_Filtering/data_f.call.gz"
    output: "snps.txt"
    message: "Extracting marker information from Lep-Map3 data file {input}"
    shell: "scripts/extract_markers.sh {input}"


rule generate_intervals:
  input:
    markers = "snps.txt",
    intervals = expand("7_Intervals/ordered.{x}.intervals", x = range(1, lg + 1))
  output: 
    intervals = "10_Anchoring/lepmap3_intervals.la"
  message: "Combining {params} Lep-Map3 interval files into single LepAnchor input {output}"
  params:
    lg = lg
  shell: 
    """
    for i in $(seq 1 {params.lg}); do
      awk -vn=$i '(NR==FNR){{map[NR-1]=$0}}(NR!=FNR){{$1=map[$1] "\t" n;print}}' {input.markers} 7_Intervals/ordered.$i.intervals >> {output.intervals}
    done
    """

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
    echo "Detected $(wc -l {output}) potentially chimeric contigs"
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
    cat {output.lift} | sort -V -k 1,1 -k 2,2n > {output.sortedlift}
    """

rule cleanmap:
  input: "10_Anchoring/liftover.sorted.la"
  output: "10_Anchoring/map_all.clean"
  log: "10_Anchoring/cleamap.log"
  message: "Running CleanMap"
  shell: "java -cp LA CleanMap map={input} > {output} 2> {log}"

rule map2bed:
  input: 
    cleanmap = "10_Anchoring/map_all.clean",
    lengths = "10_Anchoring/contigs.length",
  output: "10_Anchoring/map.bed"
  log: "10_Anchoring/map2bed.log"
  message: "Running Map2Bed"
  shell: "java -cp LA Map2Bed map={input.cleanmap} contigLength={input.lengths} > {output} 2> {log}"

rule ungrouped:
  input:
    lengths = "10_Anchoring/contigs.length",
    haplos = "10_Anchoring/fullHaplotypes50.txt",
    bedfile = "10_Anchoring/map.bed"
  output:
    bedfile = "10_Anchoring/map_extra.bed"
  message: "Finding contigs not put into chromosomes"
  params:
    chrom = lg
  shell:
    """
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {input.haplos}; cut -f 1 {input.bedfile}) > 10_Anchoring/unused_contigs.txt
    grep -w -F -f 10_Anchoring/unused_contigs.txt {input.lengths} | awk -vn={params.chrom} '{{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=n;++i) print s i}}' > 10_Anchoring/chr0.bed
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
  message: "Running PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp LA PlaceAndOrientContigs bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} keepEmptyIntervals=1 > {output} 2> {log}
    """

rule propogate:
  input: 
    placed = expand("10_Anchoring/orient_1/chr.{lgs}.la", lgs = lg_range),
    bedfile = "10_Anchoring/map.bed"
  output:
    propogated = "10_Anchoring/map_propogated.bed",
    tmp_prop = temp(expand("10_Anchoring/propogate/propogated.{lgs}.la", lgs = range(lg + 1))),
  message: "Propogating ...something"
  shell:
    """
    awk -f LA/scripts/propagate.awk {input.placed} > 10_Anchoring/tmp1.la
    awk -f LA/scripts/propagate.awk 10_Anchoring/tmp1.la > 10_Anchoring/tmp2.la
    i=2

    while ! cmp -s "10_Anchoring/tmp$i.la" "10_Anchoring/tmp$(( $i-1 )).la" ;do
	    awk -f LA/scripts/propagate.awk 10_Anchoring/tmp$i.la > 10_Anchoring/tmp$[$i+1].la
	    i=$[$i+1]
    done
    #create prop*.la
    awk '/^[^#]/{{++d[$1 "\t" $7+0 "\t" $8+0]; data[++line]=$0}}END{{for (i = 1; i <= line; ++i) {{$0=data[i];if (d[$1 "\t" $7+0 "\t" $8+0] == 1) fn="10_Anchoring/propogate/propogated."$5".la"; else if ($5==1) fn="10_Anchoring/propogate/propogated.0.la"; else fn=""; if (fn != "") print $0>fn}}}}' 10_Anchoring/tmp$i.la
    #awk '/^[^#]/{{++d[$1 "\t" $7+0 "\t" $8+0]; data[++line]=$0}}END{{for (i = 1; i <= line; ++i) {{$0=data[i];if (d[$1 "\t" $7+0 "\t" $8+0] == 1) fn="10_Anchoring/propogate/propogated."$5".la"; else if ($5==1) fn="10_Anchoring/propogate/propogated.0.la"; else fn=""; if (fn != "") print $0>fn}}' 10_Anchoring/tmp$i.la

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
    gunzip -fc {input.chain} | java -cp LA PlaceAndOrientContigs bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} keepEmptyIntervals=1 > {output.chrom} 2> {log.chrom}
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
      awk -f LA/scripts/prune.awk 10_Anchoring/orient_2/ichr.$i.la > 10_Anchoring/orient_2/ichr.${{i}}.pruned.la
    done 2> {output.pruned}
    awk -f LA/scripts/removeOverlaps.awk {input.bedfile} 10_Anchoring/orient_2/ichr.*.pruned.la > {output.cleaned}
    """

rule construct_agp:
  input:
    cleaned = "10_Anchoring/overlaps_rm.la"
  output:
    agp = "10_Anchoring/agp/chr.{lg_range}.agp",
    scaff_agp = "10_Anchoring/agp_scaffolds/chr.{lg_range}.scaffolds.agp"
  message: "Creating AGP files for linkage group {params.chrom}"
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
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {input.haplos};awk '($5!="U"){{print $6}}' {input.agp}) > {output.txt}
    grep -F -w -f {output.txt} {input.lengths} | awk '{{print $1,1,$2,1,"W",$1,1,$2,"+"}}' > {output.agp}
    cat {input.agp} {output.agp} > {output.final_agp}
    cat {input.scaff_agp} {output.agp} > {output.scaff_agp}
    """


rule build_contig_fasta:
  input:
    assembly = geno,
    agp = "10_Anchoring/REF_LA.agp"
  output:
    fasta = "10_Anchoring/Anchored.contigs.fa.gz",
  message: "Constructing final contig fasta file {input.agp}"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f LA/scripts/makefasta.awk - {input.agp} | gzip > {output.fasta}
    """

rule build_scaffold_fasta:
  input:
    assembly = geno,
    scaff_agp = "10_Anchoring/REF_LA_scaffolds.agp"
  output:
    scaff = "10_Anchoring/Anchored.scaffolds.fa.gz"
  message: "Constructing final scaffold fasta file {input.scaff_agp}"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f LA/scripts/makefasta.awk - {input.scaff_agp} | gzip > {output.scaff}
    """

rule mareymap_data:
  input:
    lift = "10_Anchoring/liftover.la",
    agp = expand("10_Anchoring/agp/chr.{lgs}.agp", lgs = lg_range)
  output: 
    mareydata = "11_MareyMaps/marey.data.gz"
  log: "11_MareyMaps/missing_scaffolds.txt"
  message: "Creating Marey map interval data"
  params:
    chrom = lg
  shell:
    """
    for c in $(seq 1 {params.chrom})
    do
      awk -vn=$c '($3==n)' {input.lift} | awk -f LA/scripts/liftover.awk 10_Anchoring/agp/chr.$c.agp - 2>> {log} | awk -vm=1 '(/LG/ && NF>=4){{if (NF==4) $5=$4;print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$5}}' | gzip
    done > {output.mareydata}
    """

rule mareymaps:
  input:
    data = "11_MareyMaps/marey.data.gz",
    agp = expand("10_Anchoring/agp/chr.{lgs}.agp", lgs = lg_range)
  output: expand("11_MareyMaps/LG.{lgs}.mareymap.png", lgs = lg_range)
  message: "Creating Marey Maps"
  shell: "Rscript LA/scripts/plot_marey.R {input.data} 10_Anchoring/agp"
