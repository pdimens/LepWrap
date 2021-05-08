rule separate_chromosomes:
    input:
        "2_Filtering/data_f.call.gz"
    output:
        "3_SeparateChromosomes/map.{lod_range}"
    log:
        "3_SeparateChromosomes/logs/map.{lod_range}.log"
    message: "Creating map for lodLimit={params.lod} >> 3_SeparateChromosomes/map.{params.lod}"
    threads: sepchrom_threads
    params:
        lod = "{lod_range}",
        dist_lod = "distortionLod=1",
    shell:
        """
        zcat {input} | java -cp LM3 SeparateChromosomes2 data=- sizeLimit=5 {informative} lodLimit={params.lod} {params.dist_lod} numThreads={threads} > {output} 2> {log}
        """

rule map_summary:
    input: expand("3_SeparateChromosomes/map.{LOD}", LOD = lod_range)
    output: "3_SeparateChromosomes/all.maps.summary"
    message: "Summarizing SeperateChromosomes2 maps >> 3_SeparateChromosomes/all.maps.summary"
    shell: "scripts/map_summary.r 3_SeparateChromosomes"

rule join_singles:
    input:
        datacall = "2_Filtering/data_f.call.gz",
        map_summ = "3_SeparateChromosomes/all.maps.summary"
    output:
        "map.master"
    log: "3_SeparateChromosomes/chosen.map"
    threads: sepchrom_threads
    message: "Joining singles"
    params:
        lod_limit = lod_lim,
        lod_diff = lod_diff,
        iterate = "iterate=1",
    shell:
        """
        echo -n -e '\nWhich map would you like to use (e.g. map.15)? map.'
        read -r
        echo "# the map chosen to use with OrderMarkers2" > {log}
        echo "map.$REPLY" >> {log}
        zcat {input.datacall} | java -cp LM3 JoinSingles2All map=3_SeparateChromosomes/map.$REPLY data=- {params.lod_limit} {params.lod_diff} {params.iterate} numThreads={threads} > {output}
        echo "Your chosen map can be found in the working directory as {output}"
        echo "A record of your choice can be found in {log}"
        sleep 5s
        """
