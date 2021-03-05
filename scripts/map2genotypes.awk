#converts phased data to "genotypes"
#usage: 
#java ... OrderMarkers2 ... outputPhasedData=1 > order_with_phase_LM3.txt
#awk [-vchr=X] [-vfullData=1] -f map2genotypes.awk order_with_phase_LM3.txt
#output columns marker name, chr, male postion, female postion, genotypes coded as "1 1", "1 2", "2 2" and 0 as missing
#providing fullData ouputs parents and pedigree...
BEGIN{
	map["00"]="1 1"
	map["01"]="1 2"
	map["10"]="2 1"
	map["11"]="2 2"
	map["0-"]="1 0"
	map["-0"]="0 1"
	map["-1"]="0 2"
	map["1-"]="2 0"
        map["--"]="0 0"
	if (chr == "")
		chr = 0
}
(/^[^#]/){
	if (!notFirst && fullData){
		notFirst = 1
		s1 =  "MARKER\tCHR\tMALE_POS\tFEMALE_POS"
		s2 =  "MARKER\tCHR\tMALE_POS\tFEMALE_POS"
		s3 =  "MARKER\tCHR\tMALE_POS\tFEMALE_POS"
		s4 =  "MARKER\tCHR\tMALE_POS\tFEMALE_POS"
		s5 =  "MARKER\tCHR\tMALE_POS\tFEMALE_POS"
		s6 =  "MARKER\tCHR\tMALE_POS\tFEMALE_POS"
		for (i = 7; i<=NF; i+=3) {
			n = length($i) / 2
			p1 = "P" (++numParents)
			p2 = "P" (++numParents)
			s1 = s1 "\t" p1 "x" p2 "\t" p1 "x" p2
			s2 = s2 "\t" p1 "\t" p2
			s3 = s3 "\t" 0 "\t" 0
			s4 = s4 "\t" 0 "\t" 0
			s5 = s5 "\t" 1 "\t" 2
			s6 = s6 "\t" 0 "\t" 0
			for (j = 1; j <= n; ++j) {
				s1 = s1 "\t" p1 "x" p2
				s2 = s2 "\tC" (++numOffspring)
				s3 = s3 "\t" p1
				s4 = s4 "\t" p2
				s5 = s5 "\t0" 
				s6 = s6 "\t0" 
			}
		}
		print s1
		print s2
		print s3
		print s4
		print s5
		print s6
	}
	s = $1 "\t" chr "\t" $2 "\t" $3
	for (i = 7; i<=NF; i+=3) {
		if (fullData) #parental data
			s = s "\t1 2\t1 2"
		n = length($i) / 2
		p1 = substr($i,1,n)
		p2 = substr($i,n+1)
		for (j = 1; j <= n; ++j)
			s = s "\t" map[substr(p1, j, 1) substr(p2, j, 1)]		
	}
	print s
}
