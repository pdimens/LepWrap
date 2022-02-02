#cuts bed based on extra cut sites
#usage: awk -f cutBed.awk map.bed cut_sites.txt >cut.bed
#       awk -f cutBed_fixN.awk contigs.length cut.bed >cut_fix.bed
#map.bed is from Map2Bed
#cut_sites is contig, start [,stop]
#both input files with non-overlapping intervals and sorted by start
BEGIN{
#        FS="\t"
        OFS="\t"
}
#store bed
(NR==FNR  && /^[^#]/){
	end[$1] = substr($3, index($3, "-") + 1) + 0
	if (!($1 in c))
		start[$1] = $2 + 0
	++c[$1]
	bed[$1,c[$1]]=$0
	beds[$1,c[$1]] = substr($2, index($2, "-") + 1) + 0
	bede[$1,c[$1]] = $3 + 0
	bedc[$1,c[$1]] = $5
	contigs[$1]
	if ($5 > maxChr)
		maxChr = $5
}
#store cutsites
(NR!=FNR){
	if (NF < 2)
		next
	#gap can be in ($3-$2+2) positions, start-1...end are the possible positions
	if (NF >= 3)
		$2 = $2 - 1
	else
		$3 = $2
	
	#remove cut if it intersects several contig boundaries
	#trim if only one...
	trim = 0
	for (i = 1; i <= c[$1]; ++i) {
		if ($2 < beds[$1,i] && $3 >= beds[$1,i]) {
			if (trim > 0) {
				print "cutBed: " $0 " removed" >"/dev/stderr"
				next
			}
			trim = i
		}
		if ($2 <= bede[$1,i] && $3 > bede[$1,i]) {
			if (trim > 0) {
				print "cutBed: " $0 " removed" >"/dev/stderr"
				next
			}
			trim = i	
		}
	}
	if (trim > 0)
		for (i = 1; i <= c[$1]; ++i) {
			if ($2 < beds[$1,i] && $3 >= beds[$1,i]) {
				printf "cutBed: " $0 " trimmed => " >"/dev/stderr"
				$3 = beds[$1,i] - 1
				print $0 >"/dev/stderr"
				
			}
			if ($2 <= bede[$1,i] && $3 > bede[$1,i]) {
				printf "cutBed: " $0 " trimmed => " >"/dev/stderr"
				$2 = bede[$1,i] + 1
				print $0 >"/dev/stderr"
			}
		}
	cut1[$1,++d[$1]] = $2
	cut2[$1,  d[$1]] = $3
}
END{
	for (contig in contigs) {
		if (d[contig] + 0 == 0) { #no cuts
			for (i = 1; i <= c[contig]; ++i)
				print bed[contig, i]
		} else { #at least one cut for contig
			j = 1
			i = 1
			st = start[contig]

			#before bed i, split zero chr
			printed = 0
			while (j <= d[contig] && cut2[contig, j] < beds[contig, i]) {
				c1 = cut1[contig, j]
				c2 = cut2[contig, j]
				for (chr = 1; chr <= maxChr; ++chr) ## put to all chromosomes
					print contig, st, c1"-"c2, "?", chr
				st = (c1+1)"-"(c2+1)
				++j
				printed = 1
			}
			if (printed == 0)
				st = st "-" beds[contig, i]

			while (i <= c[contig]) {
				#split bed i, chr in bedc				
				#printed = 0
				while (j <= d[contig] && cut1[contig, j] >= beds[contig, i] && cut2[contig, j] <= bede[contig, i]) {
					c1 = cut1[contig, j]
					c2 = cut2[contig, j]
					print contig, st, c1"-"c2, "?", bedc[contig, i]
					st = (c1+1)"-"(c2+1)
					++j
					#printed = 1
				}
				
				#after bed i, chr in bedc for first
				printed = 0
				while (j <= d[contig] && cut1[contig, j] > bede[contig, i] && (i == c[contig] || cut2[contig, j] < beds[contig, i + 1])) {
					c1 = cut1[contig, j]
					c2 = cut2[contig, j]
					if (printed == 0) {
						print contig, st, c1"-"c2, "?", bedc[contig, i]
						printed = 1
					} else
						for (chr = 1; chr <= maxChr; ++chr) ## put to all chromosomes
							print contig, st, c1"-"c2, "?", chr
					st = (c1+1)"-"(c2+1)
					++j
				}
				if (printed == 0) {
					if (i == c[contig])
						e2 = end[contig]"*"
					else
						e2 = beds[contig, i + 1] - 1
					e1 = bede[contig, i]
					print contig, st, e1"-"e2, "?", bedc[contig, i]
					st = (e1+1)"-"(e2+1)
				} else if (i == c[contig]) {
					for (chr = 1; chr <= maxChr; ++chr) ## put to all chromosomes
						print contig, st, end[contig]"*", "?", chr
	
				}

				++i
			}
		}
	}
	for (contig in d) 
		if (!(contig in contigs)){
			st = 1
			for (j = 1; j <= d[contig]; ++j) {
				c1 = cut1[contig, j]
				c2 = cut2[contig, j]
				for (chr = 1; chr <= maxChr; ++chr) ## put to all chromosomes
					print contig, st, c1"-"c2,"?",chr
				st = (c1+1)"-"(c2+1)
			}
			for (chr = 1; chr <= maxChr; ++chr) ## put to all chromosomes
				print contig, st, "N*","?",chr
		}	
}
