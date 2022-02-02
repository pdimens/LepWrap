#awk -f joinIntervals.awk contig.lengths intervals
(NR==FNR){
	l[$1]=$2
}
(NR!=FNR){
	if (!($1 in s)) {
		s[$1] = ++numContigs
		s2[numContigs] = $1
	}

	intervals1[$1, ++count[$1]] = $2
	intervals2[$1, count[$1]] = $3
	chr[$1, count[$1]] = $4
}
END{
	for (contig = 1; contig <= numContigs; ++contig) {
		c = s2[contig]
		intervals1[c, 1] = 1
		intervals2[c, count[c]] = l[c]
		for (j = 1; j < count[c]; ++j) {
			a = int(0.5 * (intervals2[c, j] + intervals1[c, j + 1]))
			intervals2[c, j] = a
			intervals1[c, j + 1] = a + 1
		}
		for (j = 1; j <= count[c]; ++j)
			print c "\t" intervals1[c, j] "\t" intervals2[c, j] "\t" chr[c, j]

	}
}
