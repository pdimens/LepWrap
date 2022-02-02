rule prep_geno:
  input: geno
  output: "8_Repeatmask/inputgenome/lepanchorinput.fa"
  message: "Preparing genome for repeat masking"
  shell:
    """
    mkdir -p 8_Repeatmask/inputgenome
    if (file {input} | grep -q compressed); then
      echo "- Assembly is compressed, creating decompressed copy: {output}"
      gunzip --stdout {input} > {output}
    else
      echo "Symlinking {input} to {output}"
      ln -srf {input} {output}
    fi
    """


rule repeatmask:
  input: "8_Repeatmask/inputgenome/lepanchorinput.fa"
  output: "8_Repeatmask/repeatmasked.fa.gz"
  log: "8_Repeatmask/Red.log"
  message: "Using Red to repeat-mask {input}"
  threads: 30
  shell:
    """
    echo "- Running Red"
    software/LepAnchor/deps/Red -gnm 8_Repeatmask/inputgenome -msk 8_Repeatmask -sco 8_Repeatmask -cnd 8_Repeatmask -rpt 8_Repeatmask > {log} 2>> {log}
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
  shell:
    """
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
  shell:
    """   
    cd 9_Chain
    ../software/LepAnchor/deps/step2.HM2 repeatmasked {threads} && rm -r repeatmasked.repeatmaskedx.result/raw.axt
    ln -sr ../{output.original} ../{output.slink}
    """
