rule clean:
    shell:
    """
    rm -rf data.call.gz data_f.call.gz maps.splitchrom map.master ordermarkers reordermarkers intervals distances distances_sexaverage 
    """

rule clean_data:
    shell:
    """
    rm -f data.call.gz data_f.call.gz
    """

rule clean_map:
    shell:
    """
    rm -rf map.master maps.splitchrom
    """

rule clean_order:
    shell:
    """
    rm -rf ordermarkers
    """

rule clean_trim:
    shell:
    """
    rm -rf ordermarkers/best.trim ordermarkers/logs/trimming
    """

rule clean_reorder:
    shell:
    """
    rm -rf reordermarkers
    """

rule clean_distance:
    shell:
    """
    rm -rf intervals distances distances_sexaverage
    """