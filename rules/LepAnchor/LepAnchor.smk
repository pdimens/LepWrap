from os import path
import glob

# PlaceAndOrientContigs #
lg = config["LepAnchor"]["lg_count"]
geno = config["LepAnchor"]["assembly"]
paf = config["LepAnchor"]["PAF_file"]
data_type = config["LepAnchor"]["PlaceAndOrientContigs"]["lepanchor_input"]
proximity = config["LepAnchor"]["proximity_file"]
haplo_limit = config["LepAnchor"]["PlaceAndOrientContigs"]["haplotype_limit"]
place_orient_extra = config["LepAnchor"]["PlaceAndOrientContigs"]["extra_params"]
# CleanMap #
cleanmap_extra = config["LepAnchor"]["CleanMap"]["extra_params"]
# Map2Bed #
map2bed_extra = config["LepAnchor"]["Map2Bed"]["extra_params"]
# Edge Trimming #
edgelen = config["LepAnchor"]["EdgeTrimming"]["edge_length"]
trimdist = config["LepAnchor"]["EdgeTrimming"]["cutoff"]
# Other #
lg_range = list(range(1,lg+1))

include: "generate_inputs.smk"
include: "mask_and_chain.smk"
include: "place_orient1.smk"
include: "place_orient2.smk"
include: "place_orient3.smk"
#include: "place_orient4.smk" 
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
