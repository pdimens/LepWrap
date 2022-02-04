rule place_orient4:
    input:
        chain = "9_Chain/chainfile.gz",
        bedfile = "10_PlaceAndOrientContigs/map.propogated2.nohaplo.bed",
        paf = paf,
        prox = proximity,
        lift = "10_PlaceAndOrientContigs/liftover.la",
        chrom = "10_PlaceAndOrientContigs/orient_1/chr.{lg_range}.la",
        chromlast = "10_PlaceAndOrientContigs/orient_3/chr.{lg_range}.la",
        liftover = expand("10_PlaceAndOrientContigs/liftover/chr.{lgs}.liftover", lgs = lg_range)
    output:
        chrom = "10_PlaceAndOrientContigs/orient_4/chr.{lg_range}.la",
        err = "10_PlaceAndOrientContigs/orient_4/errors/chr.{lg_range}.errors"
    message: "Running 4th round of PlaceAndOrientContigs for linkage group {params.chrom}"
    params:
        chrom = "{lg_range}",
        extras = place_orient_extra,
        datatype = data_type,
    threads: 5
    shell:
        """
        gunzip -fc {input.chain} | java -cp software/LepAnchor PlaceAndOrientContigs numThreads={threads} $(awk -f software/LepAnchor/scripts/pickorientation.awk {input.chrom}) bed={input.bedfile} chromosome={params.chrom} map=$MAPL chain=- paf={input.paf} proximity={input.prox} {params.datatype} {params.extras} evaluateAnchoring={input.chromlast} improveAnchoring=1 > {output.chrom} 2> {output.err}
        """

rule prune_contigblocks:
    input: "10_PlaceAndOrientContigs/orient_4/chr.{lg_range}.la"
    output: 
        chrom = "10_PlaceAndOrientContigs/pruned/chr.{lg_range}.pruned.la",
        err = temp("10_PlaceAndOrientContigs/pruned/err/chr.{lg_range}.pruned.err")
    message: "Pruning contig blocks without map support and removing overlaps"
    params:
        chrom = lg
    shell: "awk -f software/LepAnchor/scripts/prune.awk {input} > {output.chrom} 2> {output.err}"

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
        awk -f software/LepAnchor/scripts/removeOverlaps.awk {input.bedfile} {input.prunedchrom} > {output.overlaps}
        """