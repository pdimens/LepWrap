#samtools view -f 0x40 -q 20 5148.bam|awk -f umi.awk|awk -vbin=10000 -f prox10x.awk
#low memory version
BEGIN{
	if (bin == "")
		bin = 10000
	ibin = 1.0 / bin
	if (minD == "")
		minD = 2

}
function toBin(p) {
	return 1 + bin * int((p - 1) * ibin)
}
{
	p = toBin($2)
	i = $1 "\t" p "\t" $3
	if (!(i in used)) {
		snpUmi[$3,++numUmi[$3]] = $1"\t"p
		used[i]
	}
	if (!($1 in contigs)) {
		print "processing contig " $1 > "/dev/stderr"
		contigs[$1]=$2
	}
	else if ($2 > contigs[$1])
		contigs[$1]=$2

}
END{
	for (u in numUmi)
		for (i = 1; i <= numUmi[u]; ++i)
			umi[snpUmi[u,i], ++count[snpUmi[u,i]]] = u

	for (c in contigs) {
		for (p = 1; p <= contigs[c]; p+=bin) {
			delete dist
			s1 = c "\t" p
			for (i = 1; i <= count[s1]; ++i) {
				u = umi[s1, i]
				for (j = 1; j <= numUmi[u]; ++j) {
					s2 = snpUmi[u,j]
                                	++dist[s1"\t"s2]
				}
			}
		        for (d in dist)
				if (dist[d] >= minD)
                			print d "\t" dist[d]
		}
	}
}
