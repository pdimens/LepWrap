from os import path
import glob

configfile: "config.yaml"

geno = config["assembly"]
paf = config["PAF_file"]
proximity = config["proximity_file"]
lg = config["lg_count"]
lg_range = list(range(1,lg+1))
os_name = config["OS_info"]
map2bed_extra = config["extra_params_Map2Bed"]
cleanmap_extra = config["extra_params_CleanMap"]
place_orient_extra = config["extra_params_PlaceOrient"]
edgelen = config["LA_edge_length"]
trimdist = config["LA_trim_cutoff"]

rule all:
  input:
    fasta = "12_Fasta/Anchored.contigs.fa.gz",
    scaff = "12_Fasta/Anchored.scaffolds.fa.gz",
    fastaonly = "12_Fasta/Anchored.contigs.only.fa.gz",
    scaffonly = "12_Fasta/Anchored.scaffolds.only.fa.gz",
    mareydata = "13_MareyMapsUntrimmed/data.marey.gz",
    mareymaps = "13_MareyMapsUntrimmed/LepAnchor.mareymaps.pdf",
    trimmedmareymaps = "16_MareyMapsTrimmed/LepAnchor.mareymaps.pdf",
    trimmedmareydata= "16_MareyMapsTrimmed/data.marey.trimmed.gz",
    trimsummary = "15_Trim/LA.trim.summary.pdf"
  message: 
    """
    Lep-Anchor has finished. Good luck with the rest of your analyses!
    
    Output Files                 Location
    ====================================================
    anchored assemblies       |  12_Fasta/
    untrimmed marey maps      |  13_MareyMapsUntrimmed/
    updated linkage maps      |  14_NewIntervals/
    trimmed linkage maps      |  15_Trim/
    trimmed marey maps        |  16_MareyMapsTrimmed/
    """

rule repeatmask:
  input: geno
  output: "8_Repeatmask/repeatmasked.fa.gz"
  log: "8_Repeatmask/Red.log"
  message: "Using Red to repeat-mask {input}"
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
    software/LepAnchor/deps/Red -gnm . -msk 8_Repeatmask -sco 8_Repeatmask -cnd 8_Repeatmask -rpt 8_Repeatmask > {log} 2>> {log}
    echo "- Compressing repeat-masked genome from Red"
    gzip --stdout 8_Repeatmask/*.msk > {output} && rm 8_Repeatmask/*.msk
    """
    
  
rule chain_1:
  input: 
    geno = "8_Repeatmask/repeatmasked.fa.gz",
    ctrl = "software/LepAnchor/deps/all_lastz.ctl",
    scoremtx = "software/LepAnchor/deps/scoreMatrix.q"
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
        export PATH="$PATH:software/LepAnchor/deps/ubuntu"
    elif [ $OS == "centos5" ]
    then
        export PATH="$PATH:software/LepAnchor/deps/centOS5"
    elif [ $OS == "centos6" ]
    then
        export PATH="$PATH:software/LepAnchor/deps/centOS6"
    else
        echo "$OS is not recognized as one of Ubuntu, CentOS5, or CentOS6, defaulting to Ubuntu"
        export PATH="$PATH:software/LepAnchor/deps/ubuntu"
    fi
    ln -srf {input} 9_Chain/
    cd 9_Chain
    ../software/LepAnchor/deps/step1.HM2 repeatmasked {threads}
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
    OS=$(echo {params} | tr '[:upper:]' '[:lower:]')    
    echo "Using the $OS lastz/chainNet binaries"
    if [ $OS == "ubuntu" ]
    then
        export PATH="$PATH:software/LepAnchor/deps/ubuntu"
    elif [ $OS == "centos5" ]
    then
        export PATH="$PATH:software/LepAnchor/deps/centOS5"
    elif [ $OS == "centos6" ]
    then
        export PATH="$PATH:software/LepAnchor/deps/centOS6"
    else
        echo "$OS is not recognized as one of Ubuntu, CentOS5, or CentOS6, defaulting to Ubuntu"
        export PATH="$PATH:software/LepAnchor/deps/ubuntu"
    fi
    
    cd 9_Chain
    ../software/LepAnchor/deps/step2.HM2 repeatmasked {threads} && rm -r repeatmasked.repeatmaskedx.result/raw.axt
    ln -sr ../{output.original} ../{output.slink}
    """

rule extract_markers:
    input: "2_Filtering/data.filtered.lepmap3.gz"
    output: report("snps.txt", category = "Data")
    message: "Extracting marker information from Lep-Map3 data file {input}"
    shell: "scripts/extract_markers.sh {input}"


rule generate_intervals:
  input:
    markers = "snps.txt",
    intervals = expand("7_Intervals/ordered.{x}.intervals", x = range(1, lg + 1))
  output: 
    intervals = report("10_Anchoring/lepmap3_intervals.la", category = "Data")
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
  output: report("10_Anchoring/contigs.length", category = "Data")
  message: "Getting contig lengths"
  shell: "gunzip -fc {input} | awk -f software/LepAnchor/scripts/contigLength.awk > {output}"

rule find_haplotypes:
  input: "9_Chain/chainfile.gz"
  output: report("10_Anchoring/fullHaplotypes50.txt", category = "Logs")
  message: "Finding full haplotypes (potential chimeric contigs)"
  shell: 
    """
    gunzip -fc {input} | awk -f software/LepAnchor/scripts/findFullHaplotypes.awk > {output}
    echo "Detected $(wc -l {output}) potentially chimeric contigs"
    """

rule liftover:
  input: 
    chain = "9_Chain/chainfile.gz",
    intervals = "10_Anchoring/lepmap3_intervals.la",
    haplos = "10_Anchoring/fullHaplotypes50.txt"
  output: 
    lift = report("10_Anchoring/liftover.la", category = "Lifted Intervals"),
    sortedlift = report("10_Anchoring/liftover.sorted.la", category = "Lifted Intervals")
  message: "Running liftoverHaplotypes for the input maps"
  shell: 
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor LiftoverHaplotypes map={input.intervals} haplotypes={input.haplos} chain=- > {output.lift}
    cat {output.lift} | sort -V -k 1,1 -k 2,2n > {output.sortedlift}
    """

rule cleanmap:
  input: "10_Anchoring/liftover.sorted.la"
  output: "10_Anchoring/map_all.clean"
  log: report("10_Anchoring/cleamap.log", category = "Logs")
  message: "Running CleanMap"
  params:
    extras = cleanmap_extra
  shell: "java -cp software/LepAnchor CleanMap map={input} {params.extras} > {output} 2> {log}"

rule map2bed:
  input: 
    cleanmap = "10_Anchoring/map_all.clean",
    lengths = "10_Anchoring/contigs.length",
  output: "10_Anchoring/map.bed"
  log: report("10_Anchoring/map2bed.log", category = "Logs")
  message: "Running Map2Bed"
  params:
    extras = map2bed_extra
  shell: "java -cp software/LepAnchor Map2Bed map={input.cleanmap} contigLength={input.lengths} {params.extras} > {output} 2> {log}"

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
    chrom = report("10_Anchoring/orient_1/logs/chr.{lg_range}.la.err", category = "Anchoring I Logs")
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra
  threads: 3
  message: "Running PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs numThreads={threads} bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.extras} > {output} 2> {log}
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
    awk -f software/LepAnchor/scripts/propagate.awk {input.placed} > 10_Anchoring/tmp1.la
    awk -f software/LepAnchor/scripts/propagate.awk 10_Anchoring/tmp1.la > 10_Anchoring/tmp2.la
    i=2

    while ! cmp -s "10_Anchoring/tmp$i.la" "10_Anchoring/tmp$(( $i-1 )).la" ;do
	    awk -f software/LepAnchor/scripts/propagate.awk 10_Anchoring/tmp$i.la > 10_Anchoring/tmp$[$i+1].la
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
    chrom = report("10_Anchoring/orient_2/logs/ichr.{lg_range}.la.err", category = "Anchoring II Logs")
  params:
    chrom = "{lg_range}",
    extras = place_orient_extra
  message: "Running a second iteration of PlaceAndOrientContigs for linkage group {params.chrom}"
  shell:
    """
    gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs bed={input.bedfile} chromosome={params.chrom} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.extras} > {output.chrom} 2> {log.chrom}
    """

rule prune:
  input: 
    oriented = expand("10_Anchoring/orient_2/ichr.{lgs}.la", lgs = lg_range),
    bedfile = "10_Anchoring/map_propogated.bed"
  output: 
    pruned = report("10_Anchoring/orient_2/pruned.la", category = "Logs"),
    cleaned = report("10_Anchoring/overlaps_rm.la", category = "Logs")
  message: "Pruning contig blocks without map support and removing overlaps"
  params:
    chrom = lg
  shell:
    """
    for i in $(seq {params.chrom})
    do
      awk -f software/LepAnchor/scripts/prune.awk 10_Anchoring/orient_2/ichr.$i.la > 10_Anchoring/orient_2/ichr.${{i}}.pruned.la
    done 2> {output.pruned}
    awk -f software/LepAnchor/scripts/removeOverlaps.awk {input.bedfile} 10_Anchoring/orient_2/ichr.*.pruned.la > {output.cleaned}
    """

rule construct_agp:
  input:
    cleaned = "10_Anchoring/overlaps_rm.la"
  output:
    agp = report("11_AGP/contigs/chr.{lg_range}.agp", category = "Contig AGP Files"),
    scaff_agp = report("11_AGP/scaffolds/chr.{lg_range}.scaffolds.agp", category = "Scaffold AGP Files")
  message: "Creating AGP files for linkage group {params.chrom}"
  params:
    chrom = "{lg_range}"
  shell:
    """
    awk -vn={params.chrom} '($5==n)' {input.cleaned} | awk -vprefix="LG" -vlg={params.chrom} -f software/LepAnchor/scripts/makeagp_full2.awk - > {output.agp}
    awk -vn={params.chrom} '($5==n)' {input.cleaned} | awk -vprefix="LG" -vlg={params.chrom} -f software/LepAnchor/scripts/makeagp2.awk - > {output.scaff_agp}
    """

rule unused:
  input:
    lengths = "10_Anchoring/contigs.length",
    haplos = "10_Anchoring/fullHaplotypes50.txt",
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range),
  output: 
    txt = "11_AGP/not_used_final.txt",
    agp = "11_AGP/not_used.agp"
  message: "Finding unused contigs"
  shell:
    """
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {input.haplos};awk '($5!="U"){{print $6}}' {input.agp}) > {output.txt}
    grep -F -w -f {output.txt} {input.lengths} | awk '{{print $1,1,$2,1,"W",$1,1,$2,"+"}}' > {output.agp}
    """

rule build_final_agp:
  input:
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range),
    scaff_agp = expand("11_AGP/scaffolds/chr.{lgs}.scaffolds.agp", lgs = lg_range),
    unused = "11_AGP/not_used.agp",
  output:
    contig_agp = "11_AGP/lepanchor.contigs.only.agp",
    scaff_agp = "11_AGP/lepanchor.scaffolds.only.agp",
    contig_all_agp = "11_AGP/lepanchor.contigs.all.agp",
    scaff_all_agp = "11_AGP/lepanchor.scaffolds.all.agp"
  message: "Generating final AGP files"
  shell:
    """
    cat {input.agp} > {output.contig_agp}
    cat {input.scaff_agp} > {output.scaff_agp}
    cat {input.agp} {input.unused} > {output.contig_all_agp}
    cat {input.scaff_agp} {input.unused} > {output.scaff_all_agp}
    """

rule build_scaffold_only_fasta:
  input:
    assembly = geno,
    agp = "11_AGP/lepanchor.contigs.only.agp"
  output:
    fasta = "12_Fasta/Anchored.scaffolds.only.fa.gz",
  message: "Constructing final scaffold-only fasta file {output.fasta}"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f software/LepAnchor/scripts/makefasta.awk - {input.agp} | gzip > {output.fasta}
    """

rule build_scaffold_contig_fasta:
  input:
    assembly = geno,
    agp = "11_AGP/lepanchor.contigs.all.agp"
  output:
    fasta = "12_Fasta/Anchored.scaffolds.fa.gz",
  message: "Constructing final scaffold fasta file {output.fasta}"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f software/LepAnchor/scripts/makefasta.awk - {input.agp} | gzip > {output.fasta}
    """

rule build_contig_only_fasta:
  input:
    assembly = geno,
    scaff_agp = "11_AGP/lepanchor.scaffolds.only.agp"
  output:
    fasta = "12_Fasta/Anchored.contigs.only.fa.gz"
  message: "Constructing final contig-only fasta file {output.fasta}"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f software/LepAnchor/scripts/makefasta.awk - {input.scaff_agp} | gzip > {output.fasta}
    """

rule build_contig_fasta:
  input:
    assembly = geno,
    scaff_agp = "11_AGP/lepanchor.scaffolds.all.agp"
  output:
    fasta = "12_Fasta/Anchored.contigs.fa.gz"
  message: "Constructing final contig fasta file {output.fasta}"
  shell:
    """
    gunzip -fc {input.assembly} | awk -f software/LepAnchor/scripts/makefasta.awk - {input.scaff_agp} | gzip > {output.fasta}
    """

rule mareymap_data:
  input:
    lift = "10_Anchoring/liftover.la",
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range)
  output: 
    mareydata = "13_MareyMapsUntrimmed/data.marey.gz",
    sexavg = "13_MareyMapsUntrimmed/data.marey.sexavg.gz"
  log: report("13_MareyMapsUntrimmed/missing_scaffolds.txt", category = "Logs")
  message: 
    """
    Creating Marey map interval data:

    first points in uncertainty intervals  | {output.mareydata}
    midpoints in uncertainty intervals     | {output.sexavg}  
    """
  params:
    chrom = lg
  shell:
    """
    for c in $(seq 1 {params.chrom})
    do
      awk -vn=$c '($3==n)' {input.lift} | awk -f software/LepAnchor/scripts/liftover.awk 11_AGP/contigs/chr.$c.agp - | awk -vm=1 '(/LG/ && NF>=4){{if (NF==4) $5=$4;print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$5}}' | gzip
    done > {output.mareydata} 2> {log}
   
    for c in $(seq 1 {params.chrom})
    do
      awk -vn=$c '($3==n)' {input.lift} | awk -f software/LepAnchor/scripts/liftover.awk 11_AGP/contigs/chr.$c.agp - | awk -vm=1 '(/LG/ && NR>=4){{if (NF>4) s=0.5; else s=1;print $1"\t"$2"\t"$3"\t"m"\t"s*($4+$5)}}' | gzip
    done > {output.sexavg} 2> /dev/null
    """

rule mareymaps:
  input:
    data = "13_MareyMapsUntrimmed/data.marey.gz",
    sexavg = "13_MareyMapsUntrimmed/data.marey.sexavg.gz",
    agp = expand("11_AGP/contigs/chr.{lgs}.agp", lgs = lg_range)
  output: 
    indiv_plots = report(expand("13_MareyMapsUntrimmed/LG.{lgs}.mareymap.png", lgs = lg_range), category = "Marey Maps"),
    summary = report("13_MareyMapsUntrimmed/LepAnchor.mareymaps.pdf", category = "Marey Maps") ,
    sequential = report("13_MareyMapsUntrimmed/LepAnchor.sequentialmaps.pdf", category = "Sequential Maps"),
    SAsummary = report("13_MareyMapsUntrimmed/LepAnchor.sexavg.mareymaps.pdf", category = "Marey Maps Sex Avg"),
    SAsequential = report("13_MareyMapsUntrimmed/LepAnchor.sexavg.sequentialmaps.pdf", category = "Sequential Maps Sex Avg")
  message: "Creating Marey Maps"
  shell: 
    """
    Rscript software/LepAnchor/scripts/plot_marey.R {input.data} 11_AGP/contigs
    Rscript scripts/LASummary.r {input.data} true
    Rscript scripts/LASummarySexAvg.r {input.sexavg}
    """

rule generate_updated_intervals:
  input: "13_MareyMapsUntrimmed/data.marey.gz"
  output: "14_NewIntervals/LA.intervals.{lg_range}"
  message: "Splitting out LG {params.chrom} from {input}"
  params:
    chrom = "{lg_range}"
  shell:
    """
    zgrep "LG{params.chrom}\s" {input} > {output}
    """

rule trim_newintervals:
  input: "14_NewIntervals/LA.intervals.{lg_range}"
  output: 
    outfile = "15_Trim/LA.intervals.{lg_range}.trimmed",
    plot = "15_Trim/plots/LA.intervals.{lg_range}.trim.pdf"
  message: "Trimming edge clusters for {input}"
  params:
    edge = edgelen,
    dist = trimdist
  shell: "Rscript scripts/LATrim.r {input} {params.dist} {params.edge} 15_Trim"

rule merge_trimplots:
  input: expand("15_Trim/plots/LA.intervals.{lg}.trim.pdf", lg = lg_range)
  output: "15_Trim/LA.trim.summary.pdf"
  message: "Merging trimming plots into {output}"
  shell: "convert -density 200 {input} {output}"

rule merge_trimmedintervals:
  input: expand("15_Trim/LA.intervals.{lg}.trimmed", lg = lg_range)
  output: "16_MareyMapsTrimmed/data.marey.trimmed.gz"
  message: "Concatenating trimmed intervals to {output}"
  shell: "cat {input} | gzip -c > {output}"

rule plot_trimmedintervals:
  input: "16_MareyMapsTrimmed/data.marey.trimmed.gz"
  output: report("16_MareyMapsTrimmed/LepAnchor.mareymaps.pdf", category = "Trimmed Marey Maps")
  shell: "Rscript scripts/LASummary.r {input}"