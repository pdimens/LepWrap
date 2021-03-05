#script for marker binning...
BEGIN{
#ACxAG=AA,AC,AG,CG
	map["AA"] = "1 0 0 0 0 0 0"#00
	map["AC"] = "0 1 0 0 0 0 0"#01
	map["AG"] = "0 0 1 0 0 0 0"#10
	map["CG"] = "0 0 0 0 0 1 0"#11

	if (chr == "")
		 chr = 1
}
/^[^#]/{
	for (j = 7; j <= NF; ++j)
		if ($j ~ /#$/) {
			$j = substr($j, 1, length($j) - 1)
			oldNF = j
			break
		}
	if (oldNF == NF)
		next
	if (prev == "" && pedigree) {
		s1 = "CHR\tPOS"
		s2 = "CHR\tPOS"
		s3 = "CHR\tPOS"
		s4 = "CHR\tPOS"
		s5 = "CHR\tPOS"
		s6 = "CHR\tPOS"
		f = 1
		nt = 0
		for (j = 7; j <= oldNF; j+=3) {
			n = length($j) / 2
			s1 = s1 "\tF" f "\tF" f
			s2 = s2 "\t" (nt + 1) "\t" (nt + 2)
			s3 = s3 "\t0\t0"
			s4 = s4 "\t0\t0"
			s5 = s5 "\t1\t2"
			s6 = s6 "\t0\t0"
			for (i = 1; i <= n; ++i) {
				s1 = s1 "\tF" f
				s2 = s2 "\t" (nt + i + 2)
				s3 = s3 "\t" (nt + 1)
				s4 = s4 "\t" (nt + 2)
				s5 = s5 "\t0"
				s6 = s6 "\t0"
			}
			nt += n + 2
			++f
		}
		print s1 "\n" s2 "\n" s3 "\n" s4 "\n" s5 "\n" s6 
	}

	s = ""
	nt = 0
	for (j = 7; j <= oldNF; j+=3) {
		s = s "\t" map["AC"] "\t" map["AG"]
		n = length($j) / 2
		for (i = oldNF + nt + 1; i <= oldNF + nt + 4 *n; i+=4)
			s = s "\t" $i " " $(i+1) " " $(i+2) " 0 0 " $(i+3) " 0"
		nt += 4 * n
	}
	if (prev != s || FILENAME != prevFN)
		print $1 "\t" chr s
	prev = s
	prevFN = FILENAME
}
