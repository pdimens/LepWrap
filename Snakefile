from os import path
import glob
#FILENAMES = [os.path.basename(i) for i in glob.iglob("test/*.R1.fq.gz")]
#FILEBASENAMES = [i.replace('.R1.fq.gz', '') for i in FILENAMES]

#initial VCF file
vcf_file = [os.path.basename(i) for i in glob.glob("./*.vcf")]

# SeperateChromosomes2 params
lod_min = 20
lod_max = 30
lod_range = [i for i in range(lod_min, lod_max)]
#map_outs = ["maps.splichrom/map."+ str(i) for i in range(1,30)]

rule all:
    input:
        "maps.splitchrom/maps.summary.txt"
        #expand("maps.splitchrom/map.{LOD}", LOD = lod_range)

rule parentcall:
    input:
        vcf = expand("{vcffile}", vcffile = vcf_file),
        pedigree = "pedigree.txt"
    output:
        "data.call.gz"
    shell:
        "java -cp LM3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} removeNonInformative=1 | gzip > data.call.gz"

rule filtering:
    input:
        "data.call.gz"
    output:
        "data_f.call.gz"
    shell:
        """
        echo -e -n '\nSpecify your data tolerance (0.0001 to 0.01):  '
        read -r
        zcat data.call.gz | java -cp LM3 Filtering2 data=- dataTolerance=$REPLY | gzip > data_f.call.gz
        """

rule separatechromosomes:
    input:
        "data_f.call.gz"
    output:
        map = "maps.splitchrom/map.{lod_range}",
    threads: 8
    params:
        lod_lim = "lodLimit={lod_range}",
        dist_lod = "distortionLod=1"
    shell:
        "zcat {input} | java -cp LM3 SeparateChromosomes2 data=- {params.lod_lim} {params.dist_lod} numThreads={threads} > {output}" 

rule mapsummary:
    input:
        expand("maps.splitchrom/map.{LOD}", LOD = lod_range)
    output:
        "maps.splitchrom/maps.summary.txt"
    params:
        LOD_max = "{lod_max}"
    shell:
        "./scripts/map_summary.sh {params.LOD_max}"