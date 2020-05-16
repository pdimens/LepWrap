from os import path
import glob

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
        expand("distances/ordered.{lg}.distances", lg = lg_range),
        expand("distances_sexaverage/ordered.{lg}.sexavg", lg = lg_range),
        expand("intervals/ordered.{lg}.intervals", lg = lg_range)
    message:
        """
        LepMak3r is finished! Good luck with the rest of your analyses!
        """

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

rule order_markers:
    input:
        datacall = "data_f.call.gz",
        filt_map = "map.master"
    output:
        "ordermarkers/iterations/ordered.{lg_range}.{ITER}"
    log:
        run = "ordermarkers/logs/runs/ordered.{lg_range}.{ITER}.log",
        recomb = "ordermarkers/logs/recombination/ordered.{lg_range}.{ITER}.recombinations"
    message:
        """
        Ordering the markers with {params.dist_method} on linkage group: {params.chrom}, iteration: {params.iteration}
        """
    params:
        dist_method = "useKosambi=1",
        chrom = "{lg_range}",
        iteration = "{ITER}"
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.dist_method} chromosome={params.chrom} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule order_summary:
    input:
        expand("ordermarkers/iterations/ordered.{LG}.{iter}", LG = lg_range, iter= ITER)
    output:
        like = "ordermarkers/likelihoods.summary",
        recomb = "ordermarkers/recombination.summary"
    message:
        """
        Iteration likelihood summary >> {output.like}
        Recombination summary >> {output.recomb}
        """
    shell:
        """
        for LIKE in {input}; do 
            LG=$(echo $(basename $LIKE) | cut -d "." -f1,2)
            ITERUN=$(echo $LIKE | cut -d "." -f3)
            LIKELIHOOD=$(cat $LIKE | grep "likelihood = " | cut -d " " -f7)
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> {output.like}.tmp
        done
        sort {output.like}.tmp -k1,1V -k3,3nr > {output.like} && rm {output.like}.tmp
        Rscript scripts/RecombinationSummary.r ordermarkers
        """

rule find_bestlikelihoods:
    input:
        "ordermarkers/likelihoods.summary"
    output:
        "ordermarkers/best.likelihoods"
    message:
        """
        Identifying orders with best likelihoods for each LG >> {output}
        """
    shell:
        """
        NUMITER=$(find ordermarkers/iterations -maxdepth 1 -name "ordered.*.*" | cut -d "." -f3 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find ordermarkers/iterations -maxdepth 1 -name "ordered.*.*" | wc -l) 
        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p {input} | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "ordermarkers/iterations/$LIKELYMAP" >> {output}
        done
        """

rule trim_edges_clusters:
    input:
        "ordermarkers/best.likelihoods"
    output:
        "ordermarkers/best.trim/ordered.{lg_range}.trimmed"
    log:
        "ordermarkers/logs/trimming/ordered.{lg_range}.removed",
        "ordermarkers/logs/trimming/ordered.{lg_range}.trim.pdf"
    params:
        grep_lg = "ordermarkers/iterations/ordered.{lg_range}.",
        trim_threshold = "10"
    message:
        """
        Removing edge clusters >{params.trim_threshold}cM apart from the other markers in first+last 15% of {params.grep_lg}.
        """
    shell:
        """
        LG=$(grep -F {params.grep_lg} {input})
        Rscript scripts/LepMapp3rQA_single.r $(pwd) $LG {params.trim_threshold}
        """

rule trim_summary:
    input:
        expand("ordermarkers/best.trim/ordered.{lg}.trimmed", lg = lg_range)
    output:
        "ordermarkers/trim.summary"
    message:
        "Summarizing trim logs >> {output}"
    shell:
        """
        echo "# this is a summary of which markers were removed from which linkage group via trimming distant edge clusters" >> {output}
        echo -e "LG\trm_marker" >> {output}
        for each in ordermarkers/logs/trimming/ordered.*.removed ; do
            BASE=$(basename $each | cut -d "." -f1,2)
            sed -e "s/^/$BASE /" $each >> {output}.tmp 
        done
        sort -V {output}.tmp > {output} && rm {output}.tmp 
        """

rule reorder:
    input:
        datacall = "data_f.call.gz",
        filt_map = "map.master",
        trimlog = "ordermarkers/trim.summary",
        lg_order = "ordermarkers/best.trim/{trimfile}.trimmed"
    output:
        "reordermarkers/iterations/{trimfile}.{ITER}"
    log:
        run = "reordermarkers/logs/runs/{trimfile}.{ITER}.log",
        recomb = "reordermarkers/logs/recombination/{trimfile}.{ITER}.recombinations"
    message:
        """
        Reordering {input.lg_order}, iteration: {params.iteration}
        """
    params:
        dist_method = "useKosambi=1",
        eval_order="evaluateOrder=ordermarkers/best.trim/{trimfile}.trimmed",
        iteration = "{ITER}"
    threads: 2
    shell:
        """
        zcat {input.datacall} | java -cp LM3 OrderMarkers2 map={input.filt_map} data=- numThreads={threads} {params.eval_order} {params.dist_method} &> {log.run}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.run}.tmp > {output}
        grep "recombin" {log.run}.tmp > {log.recomb}
        awk '/#java/{{flag=1}} flag; /logL/{{flag=0}}' {log.run}.tmp > {log.run} && rm {log.run}.tmp
        """

rule reorder_summary:
    input:
        expand("reordermarkers/iterations/ordered.{LG}.{iter}", LG = lg_range, iter = ITER)
    output:
        like = "reordermarkers/likelihoods.summary",
        recomb = "reordermarkers/recombination.summary"
    message:
        """
        Iteration likelihood summary >> {output.like}
        Recombination summary >> {output.recomb}        
        """
    shell:
        """
        for LIKE in {input}; do 
            LG=$(echo $(basename $LIKE) | cut -d "." -f1,2)
            ITERUN=$(echo $LIKE | cut -d "." -f3)
            LIKELIHOOD=$(cat $LIKE | grep "likelihood = " | cut -d " " -f7)
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> {output.like}.tmp
        done
        sort {output.like}.tmp -k1,1V -k3,3nr > {output.like} && rm {output.like}.tmp
        Rscript scripts/RecombinationSummary.r reordermarkers
        """

rule find_bestlikelihoods2:
    input:
        "reordermarkers/likelihoods.summary"
    output:
        "reordermarkers/best.likelihoods"
    message:
        """
        Identifying ordered maps with best likelihoods for each LG >> {output}
        """
    shell:
        """
        NUMITER=$(find reordermarkers/iterations -maxdepth 1 -name "ordered.*.*" | cut -d "." -f3 | sort -V | uniq | tail -1)
        TOTALMAPS=$(find reordermarkers/iterations -maxdepth 1 -name "ordered.*.*" | wc -l) 
        for i in $(seq 1 $NUMITER $TOTALMAPS); do
            LIKELYMAP=$(sed -n ${{i}}p {input} | cut -f1,2 | awk '{{print $0, $1 "." $NF}}' | cut -d ' ' -f2)
            echo "reordermarkers/iterations/$LIKELYMAP" >> {output}
        done
        """

rule calculate_distances:
    input:
        data_call = "data_f.call.gz",
        lg = "reordermarkers/best.likelihoods"
    output:
        distance = "distances/ordered.{lg_range}.distances",
        sex_averaged = "distances_sexaverage/ordered.{lg_range}.sexavg",
        intervals = "intervals/ordered.{lg_range}.intervals"
    message:
        """
        Calculating sex-averaged marker distances and intervals for linkage group: {params.lg}
        """
    log:
        sex_averaged = "distances_sexaverage/logs/ordered.{lg_range}.sexavg.log",
        intervals = "intervals/logs/ordered.{lg_range}.int.log"
    params:
        dist_method = "useKosambi=1",
        lg = "{lg_range}"
    threads: 2
    shell:
        """
        LG=$(grep -F "reordermarkers/iterations/ordered.{params.lg}." {input.lg})
        cp $LG {output.distance}
        
        zcat {input.data_call} | java -cp LM3 OrderMarkers2 data=- evaluateOrder=$LG {params.dist_method} numThreads={threads} improveOrder=0 sexAveraged=1 &> {log.sex_averaged}.tmp
        sed -n '/\*\*\* LG \=/,$p' {log.sex_averaged}.tmp > {output.sex_averaged} 
        awk '/#java/{{flag=1}} flag; /*** LG =/{{flag=0}}' {log.sex_averaged}.tmp > {log.sex_averaged} && rm {log.sex_averaged}.tmp
        sed 's/LG \= 0/LG \= {params.lg}/g' {log.sex_averaged}

        zcat {input.data_call} | java -cp LM3 OrderMarkers2 data=- evaluateOrder=$LG {params.dist_method} numThreads={threads} calculateIntervals={output.intervals} > {log.intervals}
        #sed -n '/\*\*\* LG \=/,$p' {log.intervals}.tmp > {output.intervals} 
        #awk '/#java/{{flag=1}} flag; /*** LG =/{{flag=0}}' {log.intervals}.tmp > {log.intervals} && rm {log.intervals}.tmp
        #sed 's/LG \= 0/LG \= {params.lg}/g' {log.sex_averaged}
        """



#rule trimcheck:
#   input:
#       "ordermarkers/best.trim/trim.log",
#       expand("ordermarkers/best.trim/ordered.{lg}.trimmed", lg = lg_range)
#       #expand("ordermarkers/best.trim/ordered.{lg}.{iter}.trimmed", lg = lg_range, iter = ITER, allow_missing = True)
#   output:
#       "trim.done"
#   shell:
#       "touch {output}"

#rule reordercheck:
#    input:
#        expand("reordermarkers/ordered.{lg}.{iter}", lg = lg_range, iter = ITER)
#    output:
#        "reorder.done"
#    shell:
#        "touch {output}"