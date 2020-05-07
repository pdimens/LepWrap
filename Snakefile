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
        "intervals/int.done"
        #"reordermarkers/bestlikelihoods.txt"
        #expand("intervals/{trimfile}.intervals", trimfile = [i.split("/")[1].split(".txt")[0] for i in open("reordermarkers/bestlikelihoods.txt").read().splitlines()])
        
rule parentcall:
    input:
        vcf = expand("{vcffile}", vcffile = vcf_file),
        pedigree = "pedigree.txt"
    output:
        "data.call.gz"
    message:
        """
        Creating Lep-Map3 data file from VCF and pedigree files
        """
    shell:
        "java -cp LM3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} removeNonInformative=1 | gzip > data.call.gz"

rule filtering:
    input:
        "data.call.gz"
    output:
        "data_f.call.gz"
    message:
        """
        Filtering the data
        """
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
    message:
        """
        Creating maps for specified LOD range >> maps.splitchrom/map.LG
        """
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
    message:
        """
        Combining map summaries >> maps.splitchrom/maps.summary.txt
        """
    shell:
        "scripts/map_summary.sh {lod_max}"

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
        zcat {input.datacall} | java -cp LM3 JoinSingles2All map=maps.splitchrom/map.$REPLY data=- {params.lod_limit} {params.lod_diff} {params.iterate} numThreads={threads} > {output}
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
    message:
        """
        Ordering the markers on linkage groups 1:N
        This may take a while depending on the number of provided threads and requested iterations
        """
    params:
        dist_method = "useKosambi=1",
        chrom = "chromosome={lg_range}"
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.dist_method} {params.chrom} &> {log}
        grep -A 100000 \*\*\*\ LG\ \= {log} > {output}
        """

rule summarize_likelihoods:
    input:
        expand("ordermarkers/ordered.{LG}.{ITER}.txt", LG = lg_range, ITER = ITER)
    output:
        likelihoods = "ordermarkers/likelihoods.txt",
        sorted_likelihoods = "ordermarkers/likelihoods.sorted.txt"
    message:
        """
        Summarizing likelihoods from each iteration >> ordermarkers/likelihoods.txt
        Sorting iterations by likelihoods >> ordermarkers/likelihoods.sorted.txt
        """
    shell:
        """
        for LIKE in {input}; do 
            LG=$(echo $(basename $LIKE) | cut -d "." -f1,2)
            ITERUN=$(echo $LIKE | cut -d "." -f3)
            LIKELIHOOD=$(cat $LIKE | grep "likelihood = " | cut -d " " -f7)
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> {output.likelihoods}
        done
        sort {output.likelihoods} -k1,1V -k3,3nr > {output.sorted_likelihoods}
        """

rule find_bestlikelihoods:
    input:
        "ordermarkers/likelihoods.sorted.txt"
    output:
        "ordermarkers/bestlikelihoods.txt"
    message:
        """
        Identifying ordered maps with best likelihoods for each LG >> ordermarkers/bestlikelihoods.txt
        """
    shell:
        """
        LG=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | cut -d "." -f2 | sort -V | uniq)
        NUMITER=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | cut -d "." -f3 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | wc -l) 

        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p ordermarkers/likelihoods.sorted.txt | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "ordermarkers/$LIKELYMAP.txt" >> ordermarkers/bestlikelihoods.txt
        done
        """

rule trimming:
    input:
        "ordermarkers/bestlikelihoods.txt"
    output:
        expand("ordermarkers/best.trimmed/trimmed.{trimfile}", trimfile = [i.split("/")[1] for i in open("ordermarkers/bestlikelihoods.txt").read().splitlines()])
    params:
        trim_threshold = "10"
    log:
        "ordermarkers/best.trimmed/trimming.log",
        "ordermarkers/best.trimmed/bad_markers.txt",
        "ordermarkers/best.trimmed/trimming_plots.pdf"
    message:
        """
        Scanning the first and last 15% of markers in each LG and removing clusters >{params.trim_threshold}cM apart from the other markers. 
        """
    shell:
        "Rscript scripts/LepMapp3rQA.r $(pwd)/ordermarkers bestlikelihoods.txt {params.trim_threshold}"

rule reorder:
    input:
        datacall = "data_f.call.gz",
        filt_map = "map.master",
        lg_order = "ordermarkers/best.trimmed/trimmed.{trimfile}.txt"
    output:
        "reordermarkers/{trimfile}.{ITER}.txt"
    log:
        "reordermarkers/logs/{trimfile}.{ITER}.log"
    message:
        """
        Reordering the markers for each linkage group using the trimmed orders with the best likelihoods from initial ordering.
        This may take a while depending on the number of provided threads and requested iterations
        """
    params:
        dist_method = "useKosambi=1",
        eval_order="evaluateOrder=ordermarkers/best.trimmed/trimmed.{trimfile}.txt"
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.eval_order} {params.dist_method} &> {log}
        grep -A 100000 \*\*\*\ LG\ \= {log} > {output}
        """

rule summarize_likelihoods2:
    input:
        expand("reordermarkers/{reorder_file}", reorder_file = [os.path.basename(i) for i in glob.glob("reordermarkers/ordered.*.txt")])
        #expand("ordermarkers/ordered.{LG}.{ITER}.txt", LG = lg_range, ITER = ITER)
    output:
        likelihoods = "reordermarkers/likelihoods.txt",
        #sorted_likelihoods = "reordermarkers/likelihoods.sorted.txt"
    message:
        """
        Summarizing likelihoods from each iteration >> reordermarkers/likelihoods.txt
        Sorting iterations by likelihoods >> reordermarkers/likelihoods.sorted.txt
        """
    shell:
        """
        for LIKE in {input}; do
            LG=$(echo $(basename $LIKE) | cut -d "." -f1,2,3)
            ITERUN=$(echo $LIKE | cut -d "." -f4)
            LIKELIHOOD=$(cat $LIKE | grep "likelihood = " | cut -d " " -f7)
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> {output.likelihoods}
        done
        """

rule find_bestlikelihoods2:
    input:
        "reordermarkers/likelihoods.txt"
    output:
        sorted = "reordermarkers/likelihoods.sorted.txt",
        best = "reordermarkers/bestlikelihoods.txt"
    message:
        """
        Identifying ordered maps with best likelihoods for each LG >> reordermarkers/bestlikelihoods.txt
        """
    shell:
        """
        sort {input} -k1,1V -k3,3nr > {output.sorted}
        LG=$(find reordermarkers -maxdepth 1 -name "ordered.*.*.*.txt" | cut -d "." -f2 | sort -V | uniq)
        NUMITER=$(find reordermarkers -maxdepth 1 -name "ordered.*.*.*.txt" | cut -d "." -f4 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find reordermarkers -maxdepth 1 -name "ordered.*.*.*.txt" | wc -l) 

        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p {output.sorted} | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "reordermarkers/$LIKELYMAP.txt" >> {output.best}
        done
        """
rule link_best:
    input:
        "reordermarkers/bestlikelihoods.txt"
    output:
        dynamic(reordermarkers/best/ordered.{lg_iter}.txt)
    shell:
        """
        while read best; do
            ln -s $best reordermarkers/best/$(echo $best | basename)
        done < {input}
        """

#expand("reordermarkers/best/{lg}.txt", lg = [i.split("/")[1].split(".txt")[0] for i in open("reordermarkers/bestlikelihoods.txt").read().splitlines()])

rule intervals:
    input:
        best_lg = "reordermarkers/best/{best_reorder}.txt",
        datacall = "data_f.call.gz"
    output:
        intervals = "intervals/{best_reorder}.intervals",
        done = "intervals/int.done"
    message:
        """
        Calculating intervals for best reordered maps
        """
    threads: 2
    params:
        dist_method = "useKosambi=1"
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 data=- evaluateOrder={input.best_lg} numThreads={threads} {params.dist_method} calculateIntervals={output.intervals}
        touch {output.done}
        """