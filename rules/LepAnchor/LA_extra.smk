rule remove_haplos:
  input: 
    bedfile = "10_Anchoring/map_extra.bed",
    haplotypes = "10_Anchoring/suspected.haplotypes.after"
  output: "10_Anchoring/map.nohaplotypes.bed"
  message: "Creating bedfile with suspected haplotypes removed"
  shell: 
    """
    grep -w -v -f <(cut -f 2 {input.haplotypes}) {input.bedfile} > {output}
    #awk -f software/LepAnchor/scripts/removeHaplotypes.awk {input.bedfile} {input.haplotypes} > {output}"
    """