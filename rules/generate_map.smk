rule separate_chromosomes:
    input: "2_Filtering/data.filtered.lepmap3.gz"
    output: "3_SeparateChromosomes/LOD.{lod_range}"
    log: "3_SeparateChromosomes/logs/LOD.{lod_range}.log"
    message: "Clustering markers for lodLimit={params.lod} >> {output}"
    threads: 30
    params:
        lod = "{lod_range}",
        dist_lod = "distortionLod=1",
    shell:
        """
        zcat {input} | java -cp software/LepMap3 SeparateChromosomes2 data=- sizeLimit=5 {informative} lodLimit={params.lod} {params.dist_lod} numThreads={threads} > {output} 2> {log}
        """

rule map_summary:
    input: expand("3_SeparateChromosomes/LOD.{LOD}", LOD = lod_range)
    output: "3_SeparateChromosomes/all.LOD.summary"
    message: "Summarizing SeperateChromosomes2 maps >> {output}"
    shell: "scripts/MapSummary.r 3_SeparateChromosomes"

rule join_singles:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        map_summ = "3_SeparateChromosomes/all.LOD.summary"
    output: "LOD.master"
    log: "3_SeparateChromosomes/chosen.LOD"
    threads: 30
    message: "Joining singles to linkage groups"
    params:
        run_js2all = joinsingles,
        lod_limit = lod_lim,
        lod_diff = lod_diff,
        iterate = "iterate=1",
    shell:
        """
        echo -n -e '\nWhich map would you like to use (e.g. LOD.15)? LOD.'
        read -r
        echo -e "# the map chosen to use with OrderMarkers2\nLOD.$REPLY" > {log}
        echo "A record of your choice can be found in {log}"
        JS2A=$(echo {params.run_js2all} | tr '[:upper:]' '[:lower:]')
        if [ $JS2A == "true" ]; then
            zcat {input.datacall} | java -cp software/LepMap3 JoinSingles2All map=3_SeparateChromosomes/LOD.$REPLY data=- {params.lod_limit} {params.lod_diff} {params.iterate} numThreads={threads} > {output}
        else
            echo -e "\nSkipping JoinSingles2All and creating a symlink instead"
            ln -sr 3_SeparateChromosomes/LOD.$REPLY {output}
        fi
        sleep 2s
        """
