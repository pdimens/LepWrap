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
    # copy the binaries to the conda PATH
    cp -n software/LepAnchor/deps/ubuntu/* $(echo $PATH | cut -f1 -d":")
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