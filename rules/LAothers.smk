
rule chain_1:
  input: "8_Repeatmask/repeatmasked.fa.gz"
  output: "9_Chain/"
  message: "Running Lastz via HaploMerger2"
  threads: 30
  params:
    os = {os}
  shell:
    """
    OS=$(echo {params.os} | tr '[:upper:]' '[:lower:]')
    echo "Trying to use the {params.os} chainNet binaries"
    if [ $OS == "ubuntu" ]
    then
      export PATH="$PATH:LA/deps/chainNet_ubuntu"
    elif [ $OS == "centos5" ]
    then
      export PATH="$PATH:LA/deps/chainNet_centOS5"
    elif [ $OS == "centos6" ]
    then
      export PATH="$PATH:LA/deps/chainNet_centOS6"
    else
      echo "{params.os} is not recognized as one of Ubuntu, CentOS5, or CentOS6, defaulting to Ubuntu"
      export PATH="$PATH:LA/deps/chainNet_ubuntu"
    fi
    ln -sr {input} .
    deps/step1.HM2 repeatmasked {threads}
    rm repeatmasked.fa.gz
    """

rule chain_2:
  input: ""
  output: "9_Chain"
  message: "Running HaploMerger2 to generate the chain file"
  threads: 30
  shell:
    """
    OS=$(echo {params.os} | tr '[:upper:]' '[:lower:]')
    echo "Trying to use the {params.os} chainNet binaries"
    if [ $OS == "ubuntu" ]
    then
      export PATH="$PATH:LA/deps/chainNet_ubuntu"
    elif [ $OS == "centos5" ]
    then
      export PATH="$PATH:LA/deps/chainNet_centOS5"
    elif [ $OS == "centos6" ]
    then
      export PATH="$PATH:LA/deps/chainNet_centOS6"
    else
      echo "{params.os} is not recognized as one of Ubuntu, CentOS5, or CentOS6, defaulting to Ubuntu"
      export PATH="$PATH:LA/deps/chainNet_ubuntu"
    fi
    deps/step2.HM2 {input} {threads}
    """

rule extract_markers:
  input: "2_Filtering/data_f.call.gz"
  output: "snps.txt"
  message: "Extracting marker information from Lep-Map3 data file {input}"
  shell: "scripts/extract_markers.sh {input}"

rule generate_intervals:
  input:
    markers = "snps.txt",
    intervals = expand("7_Intervals/ordered.{x}.intervals", x = range(1,exp_lg + 1))
  output: 
    intervals = "8_Anchoring/lepmap3_intervals.la"
  message: "Combining Lep-Map3 intervals files into LepAnchor input"
  params:
    lg = {exp_lg}
  shell: 
    """
    for i in $(seq 1 {params.lg}); do
      awk -vn=$i '(NR==FNR){{map[NR-1]=$0}}(NR!=FNR){{$1=map[$1] "\t" n;print}}' {input.markers} 7_Intervals/ordered.$i.intervals >> {output.intervals}
    done
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