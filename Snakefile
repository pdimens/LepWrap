from os import path
import glob

configfile: "config.yaml"

#initial VCF file
#vcf_file = [os.path.basename(i) for i in glob.glob("./*.vcf")]

# SeperateChromosomes2 params

lod_range = list(range(config["lod_min"], config["lod_max"]+1))
lg_range = list(range(1, config["exp_lg"]+1))
ITER = list(range(1,config["iter"]+1))

include: "rules/data_prep.smk"
include: "rules/generate_map.smk"
include: "rules/order.smk"
include: "rules/trim.smk"
include: "rules/reorder.smk"
include: "rules/distances.smk"

rule all:
    input:
        expand("distances/ordered.{lg}.distances", lg = lg_range),
        expand("distances_sexaverage/ordered.{lg}.sexavg", lg = lg_range),
        expand("intervals/ordered.{lg}.intervals", lg = lg_range),
        "ordermarkers/trim.summary"
    message:
        """
        LepMak3r is finished! Good luck with the rest of your analyses!
        """