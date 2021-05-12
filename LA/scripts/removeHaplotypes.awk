#awk removeHaplotypes map.bed new_haplotypes.txt >map_nohap.bed
#removes contig regions based on haplotypes
#haplotypes can be obtained by something like this:
#sort -n -r chr*.err|awk '($NF=="haplotype" && ($1>=($4-$3+1-5000)/5000) && (!(($5 SUBSEP $6 SUBSEP $7) in h))){h[$2,$3,$4]; print}' >new_haplotypes.txt
#(haplotype score 1 per 5kb, first 5kb free)
#map.bed from Map2Bed
#chr*.err are the output from Lep-Anchor error stream, e.g. java PlaceAndOrientContigs ... chromosome=X ... 2>chrX.err
#haplotypes are put to maxChr+1
BEGIN{
	FS="\t"
	OFS="\t"
}

(NR==FNR  && /^[^#]/){
	bed[$1,++c[$1]] = $0
	if ($5 > maxChr)
		maxChr = $5
	contigs[$1]

	s1 = $2 + 0
	s2 = substr($2, index($2, "-") + 1) + 0
	if (s1 == 1)
		sm = 1
	else
		sm = int(0.5*(s1 + s2))

	e1 = $3 + 0
	e2 = substr($3, index($3, "-") + 1) + 0
	if (index($3,"*") > 0)
		em = e2
	else
		em = int(0.5*(e1 + e2))


	if (jindex[$1,sm,em] != "") {
		print "error: conflicting input bed"
		exit
	}
	jindex[$1,sm,em] = c[$1]

}
(NR!=FNR && !(($2"\t"$3"\t"$4) in processed)){

	#only process first (with highest score)
	processed[$2"\t"$3"\t"$4]

	contig = $2
	start = $3
	end = $4
	astart = $8
	aend = $9

	#find correct interval
	jmax = jindex[contig, start, end]
	if (jmax == "") {
		print "error:  contig " contig " not found"
		exit
	} else {
		$0 = bed[contig, jmax]
		$5 = maxChr + 1
		$2 = astart # "1-"
		$3 = aend # substr($3, index($3, "-"))
		bed[contig, jmax] = $0
		if (jmax > 1) {
			$0 = bed[contig, jmax - 1]
			e2 = substr($3, index($3, "-") + 1) + 0
			$3 = min(astart - 1, e2)
			bed[contig, jmax - 1] = $0
		}
		if (jmax < c[contig]) {
			$0 = bed[contig, jmax + 1]
			s1 = $2 + 0
			$2 = max(aend + 1, s1)
			bed[contig, jmax + 1] = $0
		}
	}
}
END{
	if (jmax != "")
		for (i in contigs) {
			for (j = 1; j <= c[i]; ++j)
				print bed[i, j]
		}
}
function max(a,b) {
	if (a+0 > b+0)
		return a
	return b
}
function min(a,b) {
	if (a+0 <= b+0)
		return a
	return b
}
