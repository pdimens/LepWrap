#paste snps.txt map.txt|awk -vnumSNPs=8 -vnumBps=4000 -f cleanmap.awk 
#snps.txt and map.txt should be comment free...
BEGIN{
	if (numSNPs == "")
		numSNPs = 2
	if (numBps == "")
		numBps = 10000
	OFS="\t"
}
($3 > 0 && /^[^#]/){
	if ($1 != prevScaf || $3 != prevChr) {
		if (prevScaf != "") {
			++intervals
			intervalScaf[intervals] = prevScaf
			intervalStart[intervals] = prevStart
			intervalEnd[intervals] = prevPos
			intervalChr[intervals] = prevChr
			intervalSize[intervals] = n
		}
		prevStart = $2
		n = 1
	} else
		++n
	if ($1 != prevScaf)
		pruneIntervals()

	prevScaf = $1
	prevPos = $2
	prevChr = $3

}
END{
	if (prevScaf != "") {
		++intervals
		intervalScaf[intervals] = prevScaf
		intervalStart[intervals] = prevStart
		intervalEnd[intervals] = prevPos
		intervalChr[intervals] = prevChr
		intervalSize[intervals] = n
#		print intervalScaf[intervals], intervalStart[intervals], intervalEnd[intervals], intervalChr[intervals]
	}
	pruneIntervals()
}

function pruneIntervals(      i,j,l,ii,minL)
{
	if (intervals == 0)
		return
#	print "pruning " intervals " intervals"
	if (intervals == 1) { # only one interval
		print intervalScaf[intervals], intervalStart[intervals], intervalEnd[intervals], intervalChr[intervals]
	} else {
		while (1) {
			minL = 1e100
			ii = 0
			for (i = 1; i <= intervals; ++i) {
				l = intervalEnd[i] - intervalStart[i] + 1
				if (l < minL && (l < numBps || intervalSize[i] < numSNPs)) {
					minL = l
					ii = i
				}
			}
			if (ii > 0) { # remove interval ii
#				print "delete " ii " " minL " " intervalSize[ii] " " intervalChr[ii]
				j = ii + 1
				if (ii > 1 && ii < intervals && intervalChr[ii - 1] == intervalChr[ii + 1]) { # join two neigbouring intervals
					intervalEnd[ii - 1] = intervalEnd[ii + 1]
					intervalSize[ii - 1] += intervalSize[ii + 1]
					++j
				}
				while (j <= intervals) {
					intervalScaf[ii] = intervalScaf[j]
					intervalStart[ii] = intervalStart[j]
					intervalEnd[ii] = intervalEnd[j]
					intervalChr[ii] = intervalChr[j]
					intervalSize[ii] = intervalSize[j]
					++ii
					++j
				}
				intervals = ii - 1
			} else
				break
		}
		for (i = 1; i <= intervals; ++i)
			print intervalScaf[i], intervalStart[i], intervalEnd[i], intervalChr[i]
	}
	intervals = 0
}
