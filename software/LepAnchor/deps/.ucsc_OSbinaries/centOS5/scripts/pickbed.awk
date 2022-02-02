#picks bed from based on propagate4
#awk -f pickbed.awk output_from_propagate4 map_extra.bed
(NR==FNR){
	contig = $1"\t"$2"\t"$3
	chr[contig, ++c[contig]]=$5
}
(NR!=FNR){
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
	
	contig = $1"\t"sm"\t"em

	if (contig in c && c[contig]>0) {
		for (i = 1; i <= c[contig]; ++i)
			print $1"\t"$2"\t"$3"\t?\t"chr[contig, i]
		c[contig] = 0
	} else {
		if (!(contig in c))
			warning[contig]
	}

}
END{
	for (i in warning)
		print  "contig " i " not included!" > "/dev/stderr"

}
