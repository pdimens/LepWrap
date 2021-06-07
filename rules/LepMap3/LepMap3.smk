from os import path
import glob

configfile: "config.yaml"

# load in parameters set in config.yaml
# data #
vcf = config["vcf"]
pedigree = config["pedigree"]
parentcall_extra = config["extra_params_ParentCall"]
# filtering #
data_tol=config["data_tol"]
filtering_extra = config["extra_params_Filtering"]
# separate chromosomes #
lod_max = str(config["lod_max"])
lod_range = list(range(config["lod_min"], lod_max+1))
informative = config["informative"]
sepchrom_extra = config["extra_params_SeparateChromosomes"]
# join singles #
joinsingles = config["run_joinsingles2all"]
lod_lim = config["lod_limit"]
lod_diff = config["lod_difference"]
js2a_extra = config["extra_params_JoinSingles"]
# ordering #
lg_count = config["exp_lg"]
lg_range = list(range(1, lg_count+1))
order_extra = config["extra_params_OrderMarkers"]
# trimming #
edge_len = str(config["edge_length"])
trim_thresh = str(config["trim_cutoff"])
# ordering II #
reorder_extra = config["extra_params_reOrderMarkers"]
# distances #
dist_method = config["distance_method"]

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
