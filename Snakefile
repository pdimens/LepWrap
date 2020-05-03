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
ITER = list(range(1,100+1))

rule all:
    input:
        "ordermarkers/likelihoods.txt"
        #expand("ordermarkers/logs/ordered.{LG}.{ITER}.log", LG = lg_range, ITER = list(range(1,100+1)))

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
        zcat {input} | java -cp LM3 Filtering2 data=- dataTolerance=$REPLY | gzip > data_f.call.gz
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
        "scripts/map_summary.sh {lod_max} "

rule joinsingles:
    input:
        datacall = "data_f.call.gz",
        map_summ = "maps.splitchrom/maps.summary.txt"
    output:
        "map.master"
    threads: 8
    params:
        lod_limit = "lodLimit=10",
        lod_diff = "lodDifference=2",
        iterate = "iterate=1",
    shell:
        """
        echo -n -e '\nWhich map would you like to use (e.g. map.15)? map.'
        read -r
        zcat {input.datacall} | java -cp LM3 JoinSingles2All map=maps.splitchrom/map.$REPLY data=- {params.lod_limit} {params.lod_diff} {params.iterate} numThreads={threads} > map.master
        echo 'Your filtered map can be found in the working directory'
        """

rule ordermarkers:
    input:
        datacall = "data_f.call.gz",
        filt_map = "map.master"
    output:
        "ordermarkers/ordered.{lg_range}.{ITER}.txt"
    log:
        "ordermarkers/logs/ordered.{lg_range}.{ITER}.log"
    params:
        dist_method = "useKosambi=1",
        chrom = "chromosome={lg_range}"
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.dist_method} {params.chrom} &> {log}
        grep -A 100000 \*\*\*\ LG\ \= {log} > {output}
        """

rule likelihoodsummary:
    input:
        "ordermarkers/ordered.{lg_range}.{ITER}.txt"
    output:
        "ordermarkers/likelihoods.txt"
    shell:
        """
        LG=$(echo {input} | cut -d "." -f1,2)
        ITERUN=$(echo {input} | cut -d "." -f3)
        LIKELIHOOD=$(head -1 {input} | tail -1 | cut -c 27-)
        echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> {output}
        """

rule bestlikelihoods:
    input:
        "ordermarkers/likelihoods.txt",
        "ordermarkers/ordered.{lg}.{iter}.txt"
    output:
        "ordermarkers/bestlikelihoods/ordered.{lg}.{iter}.txt"
    log:
        "ordermarkers/likelihoods.sorted.txt"
    shell:
        "scripts/bestlikelihood.sh"

rule trimming:
    input:
        "ordermarkers/bestlikelihoods/ordered.{lg}.{iter}"
    output:
        "ordermarkers/best.trimmed/ordered.{lg}.{iter}"
    log:
        "Trimming.log",
        "bad_markers.txt",
        "ordermarkers/best.trimmed/trimming.plots.pdf"
    params:
        trim_threshold = "10"
    shell:
        "Rscript scripts/LepMapp3rQA.r $(pwd)/ordermarkers/bestlikelihoods ordered {params.trim_threshold} > Trimming.log"