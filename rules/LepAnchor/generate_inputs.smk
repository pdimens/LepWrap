rule extract_markers:
    input: "2_Filtering/data.filtered.lepmap3.gz"
    output: report("snps.txt", category = "Data")
    message: "Extracting marker information from Lep-Map3 data file {input}"
    shell: "scripts/extract_markers.sh {input}"


rule generate_input_data:
  input:
    markers = "snps.txt",
    data = expand("7_Intervals/ordered.{x}.intervals", x = lg_range) if data_type == "noIntervals=0" else expand("7_Distances/ordered.{x}.distances", x = lg_range)
  output: 
    data = report("10_PlaceAndOrientContigs/lepanchor.input", category = "Data")
  message: "Combining {params} Lep-Map3 files into single LepAnchor input {output}"
  params:
    lg = lg,
    datatype = data_type
  shell: 
    """
    for i in $(seq 1 {params.lg}); do
      if [ {params.datatype} == "noIntervals=0" ]; then
        awk -vn=$i '(NR==FNR){{map[NR-1]=$0}}(NR!=FNR){{$1=map[$1] "\t" n;print}}' {input.markers} 7_Intervals/ordered.$i.intervals >> {output}
      else
        tail -n +3 7_Distances/ordered.$i.distances | awk -vn=$i '(NR==FNR){{map[NR-1]=$0}}(NR!=FNR){{$1=map[$1] "\t" n;print}}' {input.markers} - >> {output}
    fi
    done
    """


rule contiglengths:
  input: geno
  output: report("10_PlaceAndOrientContigs/contigs.length", category = "Data")
  message: "Getting contig lengths"
  shell: "gunzip -fc {input} | awk -f $CONDA_PREFIX/bin/contigLength.awk > {output}"


rule find_haplotypes:
  input: "9_Chain/chainfile.gz"
  output: report("10_PlaceAndOrientContigs/suspected.haplotypes.initial", category = "Logs")
  message: "Finding non-haplotype contigs not included in map.bed"
  shell: "gunzip -fc {input} | awk -f $CONDA_PREFIX/bin/findFullHaplotypes.awk > {output}"


rule liftover:
  input: 
    chain = "9_Chain/chainfile.gz",
    intervals = "10_PlaceAndOrientContigs/lepanchor.input",
    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes.initial"
  output: 
    lift = report("10_PlaceAndOrientContigs/liftover.la", category = "Lifted Intervals"),
    sortedlift = report("10_PlaceAndOrientContigs/liftover.sorted.la", category = "Lifted Intervals")
  message: "Running liftoverHaplotypes for the input maps"
  shell: 
    """
    gunzip -fc {input.chain} | java -cp $CONDA_PREFIX/bin/lepanchor LiftoverHaplotypes map={input.intervals} haplotypes={input.haplos} chain=- > {output.lift}
    cat {output.lift} | sort -V -k 1,1 -k 2,2n > {output.sortedlift}
    """


rule cleanmap:
  input: "10_PlaceAndOrientContigs/liftover.sorted.la"
  output: "10_PlaceAndOrientContigs/map_all.clean"
  log: report("10_PlaceAndOrientContigs/cleamap.log", category = "Logs")
  message: "Running CleanMap"
  params:
    extras = cleanmap_extra
  shell: "java -cp $CONDA_PREFIX/bin/lepanchor CleanMap map={input} {params.extras} > {output} 2> {log}"


rule map2bed:
  input: 
    cleanmap = "10_PlaceAndOrientContigs/map_all.clean",
    lengths = "10_PlaceAndOrientContigs/contigs.length",
  output: "10_PlaceAndOrientContigs/map.bed"
  log: report("10_PlaceAndOrientContigs/map2bed.log", category = "Logs")
  message: "Running Map2Bed"
  params:
    extras = map2bed_extra
  shell: "java -cp $CONDA_PREFIX/bin/lepanchor Map2Bed map={input.cleanmap} contigLength={input.lengths} {params.extras} > {output} 2> {log}"


rule ungrouped:
  input:
    lengths = "10_PlaceAndOrientContigs/contigs.length",
    haplos = "10_PlaceAndOrientContigs/suspected.haplotypes.initial",
    bedfile = "10_PlaceAndOrientContigs/map.bed"
  output:
    bedfile = "10_PlaceAndOrientContigs/map_extra.bed"
  message: "Finding contigs not put into chromosomes"
  params:
    chrom = lg
  shell:
    """
    cut -f 1 {input.lengths} | grep -v -w -F -f <(cut -f 2 {input.haplos}; cut -f 1 {input.bedfile}) > 10_PlaceAndOrientContigs/unused_contigs.txt
    grep -w -F -f 10_PlaceAndOrientContigs/unused_contigs.txt {input.lengths} | awk -vn={params.chrom} '{{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=n;++i) print s i}}' > 10_PlaceAndOrientContigs/chr0.bed
    cat {input.bedfile} 10_PlaceAndOrientContigs/chr0.bed > {output.bedfile}
    """
