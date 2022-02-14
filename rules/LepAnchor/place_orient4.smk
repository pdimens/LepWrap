rule place_orient4:
    input:
        chain = "9_Chain/chainfile.gz",
        bedfile = "10_PlaceAndOrientContigs/map.propogated2.nohaplo.bed",
        paf = paf,
        prox = proximity,
        lift = "10_PlaceAndOrientContigs/liftover.la",
        chrom = "10_PlaceAndOrientContigs/1_orient/chr.{lg_range}.la",
        chromlast = "10_PlaceAndOrientContigs/3_orient/chr.{lg_range}.la"
    output:
        chrom = "10_PlaceAndOrientContigs/4_orient/chr.{lg_range}.la",
        err = "10_PlaceAndOrientContigs/4_orient/errors/chr.{lg_range}.errors"
    message: "Running 4th round of PlaceAndOrientContigs for linkage group {params.chrom}"
    params:
        chrom = "{lg_range}",
        extras = place_orient_extra,
        datatype = data_type
    threads: 2
    shell:
        """
        gunzip -fc {input.chain} | java -cp $CONDA_PREFIX/bin/ PlaceAndOrientContigs chromosome={params.chrom} numThreads={threads} $(awk -f $CONDA_PREFIX/bin/pickorientation.awk {input.chrom}) bed={input.bedfile} map={input.lift} chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} evaluateAnchoring={input.chromlast} improveAnchoring=1 > {output.chrom} 2> {output.err}
        """

rule prune_contigblocks:
    input: "expand(10_PlaceAndOrientContigs/4_orient/chr.{lgs}.la", lgs = lg_range)
    output: 
        chrom = "10_PlaceAndOrientContigs/pruned/chr.{lg_range}.pruned.la",
    message: "Pruning contig blocks without map support and removing overlaps"
    params:
        chrom = lg
    shell: "awk -f $CONDA_PREFIX/bin/prune.awk {input} > {output.chrom} 2> {output.err}"

rule prune_post:
    input:
        bedfile = "10_PlaceAndOrientContigs/map.propogated2.bed",
        prunedchrom = expand("10_PlaceAndOrientContigs/pruned/chr.{lgs}.pruned.la", lgs = lg_range),
        prunederr = expand("10_PlaceAndOrientContigs/pruned/err/chr.{lgs}.pruned.err", lgs = lg_range)
    output: 
        overlaps = "10_PlaceAndOrientContigs/overlaps.removed.la",
        pruned = "10_PlaceAndOrientContigs/pruned.la"
    message: "Removing overlaps"
    threads: 1
    shell:
        """
        cat {input.prunederr} > {output.pruned}
        awk -f $CONDA_PREFIX/bin/removeOverlaps.awk {input.bedfile} {input.prunedchrom} > {output.overlaps}
        """