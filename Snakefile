from os import path
import glob
#FILENAMES = [os.path.basename(i) for i in glob.iglob("test/*.R1.fq.gz")]
#FILEBASENAMES = [i.replace('.R1.fq.gz', '') for i in FILENAMES]

#initial VCF file
vcf_file = [os.path.basename(i) for i in glob.glob("./*.vcf")]

# SeperateChromosomes2 params
lod_min = 20
lod_max = 40
lod_range = list(range(lod_min, lod_max+1))
exp_lg = 24
lg_range = list(range(1, exp_lg+1))

rule all:
    input:
        expand("ordermarkers/ordered.{LG}.{ITER}.txt", LG = lg_range, ITER = list(range(1,100+1)))

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
        map = "maps.splitchrom/map.{lod_range}"
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
    shell:
        "./scripts/map_summary.sh {lod_max} "

rule joinsingles:
    input:
        map = "maps.splitchrom/map.{lod_range}"
    output:
        "map.{lod_range}.master"
    threads: 8
    params:
        lod_limit = "lodLimit=10",
        lod_diff = "lodDifference=2",
        iterate = "iterate=1",
    shell:
        """
        echo -n -e '\nWhich map would you like to use (e.g. map.15)? map.'
        read -r
        zcat data_f.call.gz | java -cp LM3 JoinSingles2All map=maps.splitchrom/map.$REPLY data=- {params.lod_limit} {params.lod_diff} {params.iterate} numThreads={threads} > map.$REPLY.master
        echo 'Your filtered map can be found in the working directory'
        """

rule ordermarkers:
    input:
        expand("map.{LOD}.master", LOD = lod_range)
    output:
        logfile = expand("ordermarkers/logs/ordered.{LG}.{ITER}.log", LG = lg_range, ITER = list(range(1,100+1))),
        lgfile = expand("ordermarkers/ordered.{LG}.{ITER}.txt", LG = lg_range, ITER = list(range(1,100+1)))
    params:
        dist_method = "useKosambi=1",
        chrom = "chromosome={LG}"
    threads: 2
    shell:
        """
        zcat data_f.call.gz | java -cp LM3 OrderMarkers2 map={input} data=- numThreads={threads} {params.dist_method} {params.chrom} &> {output.logfile}
        grep -A 100000 \*\*\*\ LG\ \= {output.logfile} > {output.lgfile}
        """