#awk -f removeOverlaps map.bed chr*.la >no_overlap_allchr.la
#map.bed from Map2Bed
#chr*.la from PlaceAndOrientContigs
#removes overlap and trims contig cut positions
#may leave some orverlap if two alignments are overlapping...
#(assumed tab as white space separator, now FS="\t" commented)
BEGIN{
#	FS="\t"
	OFS="\t"
}

(NR==FNR  && /^[^#]/){
	bed[$1,++count[$1]] = $0
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
		print "error: conflicting input bed" >"/dev/stderr"
		print $0 >"/dev/stderr"
		exit 1
	}
	jindex[$1,sm,em] = count[$1]
}

(NR!=FNR && /^[^#]/){
	laoutput[++lines] = $0

	origs = $7 + 0 #remove alt start
	alts = substr($7, index($7, "/") + 1) + 0
	orige = $8 + 0 #remove alt end
	alte = substr($8, index($8, "/") + 1) + 0


	#find corresponding interval from bed
	jmax = jindex[$1,origs,orige]
	if (jmax == "") {
		print "error:  contig " $1 " not found" >"/dev/stderr"
		print $0 >"/dev/stderr"
#		exit 1
	}
	else {
		if ($9 != "null" && ($10 == "null" || $10 ~ /^chain/)) { #linked to prev (fixed end point)
			if (index($4, "-") > 0) #- orientation
				cute[$1,jmax] = max($3, alte)
			else
				cuts[$1,jmax] = min($2, alts)
		}
		if ($12 != "null" && ($13 == "null" || $13 ~ /^chain/)) { #linked to next (fixed end point)
			if (index($4, "-") > 0) #- orientation
				cuts[$1,jmax] = min($2, alts)
			else
				cute[$1,jmax] = max($3, alte)
		}
	}
}
END{
#	if (jmax == "")
#		exit

	#propagate cut points
	for (contig in contigs) {
		if (count[contig] > 1) { #if only one region, no need to update...
			for (j=1;j<=count[contig]; ++j) {
				if (j > 1 && cuts[contig, j] == "" && cute[contig, j - 1] != "") {
					cuts[contig, j] = cute[contig, j - 1] + 1
				}
				if (j < count[contig] && cute[contig, j] == "" && cuts[contig, j + 1] != "") {
					cute[contig, j] = cuts[contig, j + 1] - 1
				}
			}

		}
	}

	


	trimmed = 0
	#print trimmed result
	for (i = 1; i <= lines; ++i) {
		$0 = laoutput[i]

		#find corresponding interval from bed
		origs = $7 + 0 #remove alt start
		orige = $8 + 0 #remove alt end
		jmax = jindex[$1,origs,orige]

		if (jmax == "")
			print
		else {

			if ($9 == "null" || ($10 != "null" && !($10 ~ /^chain/))) { #not linked to prev (with fixed end point)
				if (index($4, "-") > 0) {#- orientation
					ce = cute[$1,jmax] 
					if (ce != "") {
						$3 = ce
						++trimmed
					}
				} else {
					cs = cuts[$1,jmax]
					if (cs != "") {
						$2 = cs
						++trimmed
					}
				}
			}
			if ($12 == "null" || ($13 != "null" && !($13 ~ /^chain/))) { #not linked to next (with fixed end point)
				if (index($4, "-") > 0) {#- orientation
					cs = cuts[$1,jmax]
					if (cs != "") {
						$2 = cs
						++trimmed
					}
				} else {
					ce = cute[$1,jmax] 
					if (ce != "") {
						$3 = ce
						++trimmed
					}
				}
			}
			print 
		}
	}
	print "trimmed " trimmed " contig ends" >"/dev/stderr"
}

function min(a,b) {
	if (a+0 <= b+0)
		return a
	return b
}
function max(a,b) {
	if (a+0 > b+0)
		return a
	return b
}
#intersection size of two intervals
function intersection(s1, e1, s2, e2)
{
	return min(e1, e2) - max(s1, s2) + 1
}
