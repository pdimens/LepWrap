from os import path
import glob

# load in parameters set in the config file
# data #
vcf = config["LepMap3"]["vcf"]
pedigree = config["LepMap3"]["pedigree"]
parentcall_extra = config["LepMap3"]["ParentCall2"]["extra_params"]
# filtering #
data_tol=config["LepMap3"]["Filtering2"]["data_tol"]
filtering_extra = config["LepMap3"]["Filtering2"]["extra_params"]
# separate chromosomes #
lod_max = str(config["LepMap3"]["SeperateChromosomes2"]["lod_max"])
lod_range = list(range(config["LepMap3"]["SeperateChromosomes2"]["lod_min"], config["LepMap3"]["SeperateChromosomes2"]["lod_max"]+1))
informative = config["LepMap3"]["SeperateChromosomes2"]["informative"]
sepchrom_extra = config["LepMap3"]["SeperateChromosomes2"]["extra_params"]
# join singles #
joinsingles = config["LepMap3"]["JoinSingles2ALL"]["run"]
lod_lim = config["LepMap3"]["JoinSingles2ALL"]["lod_limit"]
lod_diff = config["LepMap3"]["JoinSingles2ALL"]["lod_difference"]
js2a_extra = config["LepMap3"]["JoinSingles2ALL"]["extra_params"]
# ordering #
lg_count = config["LepMap3"]["OrderMarkers2"]["exp_lg"]
lg_range = list(range(1, lg_count+1))
order_extra = config["LepMap3"]["OrderMarkers2"]["extra_params"]
# trimming #
edge_len = str(config["LepMap3"]["EdgeTrimming"]["edge_length"])
trim_thresh = str(config["LepMap3"]["EdgeTrimming"]["cutoff"])
# ordering II #
reorder_extra = config["LepMap3"]["ReOrderMarkers2"]["extra_params"]
# distances #
dist_method = config["LepMap3"]["CalculateDistances"]["distance_method"]

include: "prepare_data.smk"
include: "generate_map.smk"
include: "order.smk"
include: "trim.smk"
include: "reorder.smk"
include: "distances.smk"

rule all:
    input:
        expand("7_Distances/ordered.{lg}.distances", lg = lg_range),
        expand("7_DistancesSexAverage/ordered.{lg}.sexavg", lg = lg_range),
        expand("7_Intervals/ordered.{lg}.intervals", lg = lg_range),
        "4_OrderMarkers/recombination/recombination.summary",
        "5_Trim/trim.summary",
        "6_OrderMarkers/recombination/recombination.summary",
    message:
        """
        Lep-Map3 has finished. Good luck with the rest of your analyses!
        Output Files:
        =============
        final linkage maps  | 7_Distances/ordered.*.distances
        sex-averaged maps   | 7_DistancesSexAverage/ordered.*.sexavg
        map-intervals       | 7_Intervals/ordered.*.intervals
        """
