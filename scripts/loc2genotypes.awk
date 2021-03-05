#converts a loc file (JoinMap) to genotypes, that can be converted back post file
#should handle windows end-of-line characters as well
#awk -f locsingle.awk file.loc|awk -f loc2genotypes.awk|awk -f genotypes2post.awk |java -cp Lep-MAP3/bin ... data=- ...
BEGIN{
	map2["-"] = 0
	#<lmxll>
	map2["l"] = 1
	map2["m"] = 2
	#<nnxnp>
	map2["n"] = 1
	map2["p"] = 2
	#<hkxhk>
	map2["h"] = 1
	map2["k"] = 2
	#<efxeg>
	map2["e"] = 1
	map2["f"] = 2
	map2["g"] = 3
	#<abxcd>
	map2["a"] = 1
	map2["b"] = 2
	map2["c"] = 3
	map2["d"] = 4
	for (i in map2) {
		map[i] = map2[i]
		map[toupper(i)] = map2[i]
	}
	delete map2
	markers["<lmxll>"]
	markers["<nnxnp>"]
	markers["<hkxhk>"]
	markers["<efxeg>"]
	markers["<abxcd>"]
	markers["<abxab>"]
	markers["<abxaa>"]
	markers["<aaxab>"]
}
($2 in markers) {
	if ($NF == "\r")	#windows end of line
		--NF
	++line
	start = 3
	if ($3 ~ /{.*}/ || $3 ~ /\(.*\)/)
		start = 4
#	if (start == 3 && !($3 ~ /{.*}/) {
#		print "Error: missing phase on data" >/dev/stderr
#		exit(-1)
#	}

	if (line == 1) { ## print pedigree
		s = "CHR\tPOS"
		for (i = start; i <= NF + 4; ++i)
			s = s "\tF" 
		print s
		s = "CHR\tPOS\tGP1\tGP2\tP1\tP2"
		for (i = start; i <= NF; ++i)
			s = s "\t"  (i - start + 1)
		print s
		s = "CHR\tPOS\t0\t0\tGP1\tGP2"
		for (i = start; i <= NF; ++i)
			s = s "\tP1"
		print s
		s = "CHR\tPOS\t0\t0\t0\t0"
		for (i = start; i <= NF; ++i)
			s = s "\tP2"
		print s
		s = "CHR\tPOS\t1\t1\t1\t2"
		for (i = start; i <= NF; ++i)
			s = s "\t0"
		print s
		s = "CHR\tPOS"
		for (i = start; i <= NF + 4; ++i)
			s = s "\t0" 
		print s
	}
	s = $1 "\t" line
	if (start == 4) {
		if ($3 ~ /{0/)
			s = s "\t" map[substr($2, 2, 1)] " " map[substr($2, 2, 1)]
		else  if ($3 ~ /{1/)
			s = s "\t" map[substr($2, 3, 1)] " " map[substr($2, 3, 1)]
		else 	
			s = s "\t0 0"

		if ($3 ~ /{.0/)
			s = s "\t" map[substr($2, 5, 1)] " " map[substr($2, 5, 1)]
		else  if ($3 ~ /{.1/)
			s = s "\t" map[substr($2, 6, 1)] " " map[substr($2, 6, 1)]
		else 	
			s = s "\t0 0"

	} else
		s = s "\t0 0\t0 0"
	s = s "\t" map[substr($2, 2, 1)] " " map[substr($2, 3, 1)] "\t" map[substr($2, 5, 1)] " " map[substr($2, 6, 1)]

	for (i = start; i <= NF; ++i)
		s = s "\t" map[substr($i, 1, 1)] " " map[substr($i, 2, 1)]
	print s
}
END{
	
}
