from os import path
import glob

# PlaceAndOrientContigs #
data_type = config["lepanchor_input"]
geno = config["assembly"]
paf = config["PAF_file"]
proximity = config["proximity_file"]
haplo_limit = config["haplotype_limit"]
place_orient_extra = config["extra_params_PlaceOrient"]
# Map2Bed #
map2bed_extra = config["extra_params_Map2Bed"]
# CleanMap #
cleanmap_extra = config["extra_params_CleanMap"]
# Edge Trimming #
edgelen = config["LA_edge_length"]
trimdist = config["LA_trim_cutoff"]
# Other #
lg = config["lg_count"]
lg_range = list(range(1,lg+1))

include: "generate_inputs.smk"
include: "mask_and_chain.smk"
include: "place_orient.smk"
include: "place_orient_ii.smk"
include: "place_orient_iii.smk"
include: "build_agp.smk"
include: "build_fasta.smk"
include: "mareymaps_untrimmed.smk"
include: "trim_edges.smk"

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
