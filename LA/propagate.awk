#awk -f propagate chr[0-9]*.la
BEGIN{
	if (scoreDiff == "")
		scoreDiff = 2	#how much score is needed to assign contigs to chrs (2 pacbio reads)
	OFS="\t"
}
(/^[^#]/){
	if ($5 > maxChr)
		maxChr = $5
	data[++numLines]=$0
	contig = $1 "\t" $7+0 "\t" $8+0
	++count[contig]
}

#processing of data[line]
#Find contigs without chromosome linked to contigs in chromosomes (mode=1)
#Remove linked contigs in wrong chromosomes (mode=2)
#print contigs (mode != 1 && mode != 2)
function processLine(line, mode)
{
	$0 = data[line]
	contig = $1 "\t" $7+0 "\t" $8+0
	link1 = $9
	link2 = $12
	$0 = data[line - 1]

	if ($1 == link1)
		prevContig = $1 "\t" $7+0 "\t" $8+0
	else
		prevContig = "null"

	$0 = data[line + 1]
	if ($1 == link2)
		nextContig = $1 "\t" $7+0 "\t" $8+0
	else
		nextContig = "null"

	$0 = data[line]

	if (mode == 1) {
		if (count[contig] > 1) {
			if (count[prevContig] == 1 && $11 + $15 >= score[contig,$5])
				score[contig,$5]=$11 + $15
			if (count[nextContig] == 1 && $14 + $15 >= score[contig,$5])
				score[contig,$5]=$14 + $15
		}
	} else if (mode == 2){
	#remove extra links
		if (count[prevContig] > 1 && (prevContig in pchr) && pchr[prevContig] != $5) {
			$0 = data[line]
			$9 = "null"
			$10 = "null"
			$11 = 0
			data[line] = $0

			$0 = data[line - 1]
			$12 = "null"
			$13 = "null"
			$14 = 0
			data[line - 1] = $0
		}
		if (count[nextContig] > 1 && (nextContig in pchr) && pchr[nextContig] != $5) {
			$0 = data[line]
			$12 = "null"
			$13 = "null"
			$14 = 0
			data[line] = $0

			$0 = data[line + 1]
			$9 = "null"
			$10 = "null"
			$11 = 0
			data[line + 1] = $0
		}
	} else {
		#print ouput
		if (count[contig] == 1 || !(contig in pchr) || pchr[contig] == $5)
			print
	}
}

END{

	for (line = 1; line <= numLines; ++line) {
		processLine(line, 1)
	}

#Find unique links
#calculate propagated chromosomes, pchr
	for (nu in count) {
		if (count[nu] > 1) {
			bestScore  = 0
			bestScore2 = 0
			bestChr = 0
			for (chr = 1; chr <= maxChr; ++chr) {
				s = score[nu,chr]
				if (s >= bestScore) {
					bestScore2 = bestScore
					bestScore = s
					bestChr = chr
				} else if (s >= bestScore2)
					bestScore2 = s
			}
			if (bestScore > 0 && bestScore2 + scoreDiff <= bestScore) { # at least scoreDiff difference
#				print nu "\t?\t" bestChr
				pchr[nu] = bestChr
			} else {
				notUsed2[nu]
			}
		}
	}
	#remove extra links
	for (line = 1; line <= numLines; ++line)
		processLine(line, 2)

	#print final result
	for (line = 1; line <= numLines; ++line)
		processLine(line, 3)

}

