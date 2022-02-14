rule build_scaffold_only_fasta:
  input:
    assembly = geno,
    agp = "11_AGP/lepanchor.contigs.only.agp"
  output: "12_Fasta/Anchored.scaffolds.only.fa.gz"
  message: "Constructing final scaffold-only fasta file {output}"
  shell: "gunzip -fc {input.assembly} | awk -f $CONDA_PREFIX/bin/makefasta.awk - {input.agp} | gzip > {output}"


rule build_scaffold_contig_fasta:
  input:
    assembly = geno,
    agp = "11_AGP/lepanchor.contigs.all.agp"
  output: "12_Fasta/Anchored.scaffolds.fa.gz"
  message: "Constructing final scaffold fasta file {output}"
  shell: "gunzip -fc {input.assembly} | awk -f $CONDA_PREFIX/bin/makefasta.awk - {input.agp} | gzip > {output}"


rule build_contig_only_fasta:
  input:
    assembly = geno,
    scaff_agp = "11_AGP/lepanchor.scaffolds.only.agp"
  output: "12_Fasta/Anchored.contigs.only.fa.gz"
  message: "Constructing final contig-only fasta file {output}"
  shell: "gunzip -fc {input.assembly} | awk -f $CONDA_PREFIX/bin/makefasta.awk - {input.scaff_agp} | gzip > {output}"


rule build_contig_fasta:
  input:
    assembly = geno,
    scaff_agp = "11_AGP/lepanchor.scaffolds.all.agp"
  output: "12_Fasta/Anchored.contigs.fa.gz"
  message: "Constructing final contig fasta file {output}"
  shell: "gunzip -fc {input.assembly} | awk -f $CONDA_PREFIX/bin/makefasta.awk - {input.scaff_agp} | gzip > {output}"