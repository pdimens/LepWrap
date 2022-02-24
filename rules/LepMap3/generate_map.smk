rule separate_chromosomes:
    input: "2_Filtering/data.filtered.lepmap3.gz"
    output: "3_SeparateChromosomes/LOD.{lod_range}"
    log: "3_SeparateChromosomes/logs/LOD.{lod_range}.log"
    message: "Clustering markers for lodLimit={params.lod} >> {output}"
    threads: 40
    params:
        lod = "{lod_range}",
        extra = sepchrom_extra,
    shell:
        """
        zcat {input} | java -cp $CONDA_PREFIX/bin/lepmap3 SeparateChromosomes2 lodLimit={params.lod} data=- {params.extra} {informative} numThreads={threads} > {output} 2> {log}
        """


rule map_summary:
    input: expand("3_SeparateChromosomes/LOD.{LOD}", LOD = lod_range)
    output: "3_SeparateChromosomes/all.LOD.summary"
    message: "Summarizing SeperateChromosomes2 maps >> {output}"
    shell: "MapSummary.r 3_SeparateChromosomes"


rule choose_map:
    input: "3_SeparateChromosomes/all.LOD.summary"
    output: "3_SeparateChromosomes/chosen.LOD"
    message: "Examine {input} and decide on a map of a given LOD limit before proceeding"
    shell:
        """
        echo -n -e '\nWhich map would you like to use (e.g. LOD.15)? LOD.'
        read -r
        echo -e "# the map chosen to use with OrderMarkers2\n3_SeparateChromosomes/LOD.$REPLY" > {output}
        echo "A record of your choice can be found in {output}"
        sleep 2s      
        """


rule join_singles:
    input:
        datacall = "2_Filtering/data.filtered.lepmap3.gz",
        map_choice = "3_SeparateChromosomes/chosen.LOD"
    output: "LOD.master"
    threads: 40
    message: "Joining singles to linkage groups"
    params:
        run_js2all = joinsingles,
        lod_limit = lod_lim,
        lod_diff = lod_diff,
        extra = js2a_extra
    shell:
        """
        JS2A=$(echo {params.run_js2all} | tr '[:upper:]' '[:lower:]')
        THEMAP=$(tail -1 {input.map_choice})
        if [ $JS2A == "true" ]; then
            zcat {input.datacall} | java -cp $CONDA_PREFIX/bin/lepmap3 JoinSingles2All map=$THEMAP data=- {params.extra} {params.lod_limit} {params.lod_diff} numThreads={threads} > {output}
        else
            echo -e "\nSkipping JoinSingles2All and creating a symlink to $THEMAP instead"
            ln -sr $THEMAP {output}
        fi
        sleep 2s
        """