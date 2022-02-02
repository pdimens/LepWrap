#chains paf files
#assumes paf is sorted on first column (query string) as in output of minimap2
#awk -f chainpaf.awk aln20DP.paf 
BEGIN{
#not implemented yet
#	if (maxDistance == "")
#		maxDistance=2000

	if (gapPenalty == "")
		gapPenalty=-1
	if (gapOpen == "")
		gapOpen=-12
	OFS="\t"
	FS="\t"
}
{
	if (prev == "" || $1 == prev)
		data[++lines] = $0
	else {
		tmp = $0
		processData()
		lines=1
		data[lines] = tmp
		$0 = tmp
	}
	prev = $1
}
END{
	processData()
}
function gapScore(end, start) {
	if (start >= end) {
		#print "error"
		#exit
		return 0
	}
	return gapOpen + gapPenalty  * (end - start)
}

function max(a, b) {
	if (a >= b)
		return a
	return b
}
function min(a, b) {
	if (a < b)
		return a
	return b
}


function processData(i,j,k,n,gl1,gl2,gp1,gp2,start1,start2,end1,end2,bj,starts,cstarts,startsi,other,cother,best,prev,deleted,ret,mydata,s,cigar,score,trim)
{
	for (i = 1; i <= lines; ++i) {
		$0=data[i]
		starts[i] = $3
		startsi[$3, ++cstarts[$3]]=i
	}
	n = asort(starts)

	#sort data based on query start
	#and separate by target
	for (i = 1; i <= lines; ++i) {
		start = starts[i]
		$0 = data[startsi[start, cstarts[start]]]
		s = 0
		for (j = 13; j <= NF; ++j) {
			if ($j ~ /^AS:i:/) {
				s = substr($j, 6)+0
				break
			}
		}
		if (s <= 0) #do not store chains with score <= 0
			continue
		--cstarts[start]
		other[$6,$5]
		mydata[$6,$5,++cother[$6,$5]] = $0
	}

	for (o in other) {
		if (cother[o] == 1)
			print mydata[o,1]
		else { # find max score non-overlapping chain...
		first = 1
		while (cother[o] > 1) {
			for (j = 1; j <= cother[o]; ++j) {
				$0 = mydata[o,j]
				if (first)
					print
				for (i = 13; i <= NF; ++i) {
					if ($i ~ /^AS:i:/)
						score[j] = substr($i, 6)+0
					if ($i ~ /^cg:Z:/)
						cigar[j] = substr($i, 6)
				}
				best[j] = score[j]
				prev[j] = j
				trim[j] = 0 #trim cigar if needed
			}
			first = 0

			for (j = 1; j < cother[o]; ++j) {
				$0 = mydata[o,j]
				end1 = $4
				if ($5 == "+") {
					end2 = $9
					for (k = j + 1; k <= cother[o]; ++k) {
						$0 = mydata[o,k]
						#($5 == "+"), both alignments are in the same orientation
						start1 = $3
						start2 = $8
						#if (start1 >= end1 && start2 >= end2) { #non-overlapping, increasing
						ov = max(end1 - start1, end2 - start2) #overlap
						if (ov <= 0 || ov <= maxTrim1(cigar[k])) {
							if (ov <= 0) {
								gp1 = gapScore(start1, end1)
								gp2 = gapScore(start2, end2)
							}
							else {
								gp1 = gapScore(start1 + ov, end1)
								gp2 = gapScore(start2 + ov, end2) - score[k] * ov /  min($9 - $8, $4 - $3)
							}
							if (best[j] + gp1 <= 0)
								break
							nb = best[j] + score[k] + gp1 + gp2
							if (nb > best[k]) {
								best[k] = nb
								prev[k] = j
								trim[k] = ov
							}
						}
					}
				} else { #"-" orientation
					end2 = $8
					for (k = j + 1; k <= cother[o]; ++k) {
						$0 = mydata[o,k]
						#($5 == "-"), both alignments are in the same orientation
						start1 = $3
						start2 = $9
						#if (start1 >= end1 && start2 <= end2) { #non-overlapping, increasing (decreasing orientation -)
                                              ov = max(end1 - start1, start2 - end2) #overlap
						if (ov <= 0 || ov <= maxTrim1(cigar[k])) {
							if (ov <= 0) {
								gp1 = gapScore(start1, end1)
								gp2 = gapScore(end2, start2)
							}
							else {
                                                      	gp1 = gapScore(start1 + ov, end1)
                                                      	gp2 = gapScore(end2, start2 - ov) - score[k] * ov /  min($9 - $8, $4 - $3)
							}
							if (best[j] + gp1 <= 0)
								break
							nb = best[j] + score[k] + gp1 + gp2
							if (nb > best[k]) {
								best[k] = nb
								prev[k] = j
                                                            	trim[k] = ov
							}
						}
					}
				}
			}
			while (1) {
				bj = 1
				for (j = 2; j <= cother[o]; ++j) {
					if ((best[j] > best[bj] && !(j in deleted)) || (bj in deleted))
						bj = j
				}
				n = 1
				ret[1] = bj;
				while (prev[bj] != bj) {
					bj = prev[bj] 
					ret[++n] = bj
				}
				s = 0
				for (i = 1; i <= n; ++i)
					if (ret[i] in deleted) {
						s = 1
						break
					}
				if (s)
					break
					
				for (i = 1; i <= n; ++i)
					deleted[ret[i]]
					
				if (n == 1)
					continue

				$0 = mydata[o,ret[n]]
				start1 = $3
				if ($5 == "+") {
					start2 = $8
				} else {
					end2 =  $9
				}

				$0 = mydata[o,ret[1]]
				end1 = $4
				if ($5 == "+") {
					end2 = $9
				} else {
					start2 =  $8
				}
				s = $1 "\t" $2 "\t" start1 "\t" end1 "\t" $5 "\t" $6 "\t" $7 "\t" start2 "\t" end2 "\t0\t0\t0\tAS:i:" int(best[ret[1]]) "\t"  "CH:i:" n "\tcg:Z:"


				$0 = mydata[o,ret[n]]
				end1 = $4
				if ($5 == "+")
					end2 = $9
				else
					end2 = $8

				#construct combined cigar
				c = cigar[ret[n]]

				for (i = n - 1; i >= 1; --i) {
					j = ret[i]
					$0 = mydata[o,j]

					gl1 = $3 - end1
					end1 = $4
					if ($5 == "+") {
						gl2 = $8 - end2
						end2 = $9
					}
					else {
						gl2 = end2 - $9
						end2 = $8
					}
					if (trim[j] > 0) {
						gl1 += trim[j]
						gl2 += trim[j]
					}

					if (gl1 < 0 || gl2 < 0)  {
						print "\nerror:negative gap" gl1 gl2
						exit
					}
					if (gl1 > 0)
						c = c gl1 "I"
					if (gl2 > 0)
						c = c gl2 "D"
					if (trim[j] > 0)
						c = c "" trim1(cigar[j], trim[j])
					else	
						c = c "" cigar[j]

				}
				if (n > 1)
					print s c
			} #while (1)
				
			k = 1
			for (j = 1; j <= cother[o]; ++j) {
				mydata[o, k] = mydata[o, j]
				if (!(j in deleted))
					++k;
			}
			cother[o] = k - 1
			delete(deleted)
		}
		}
	}
}
#how much you can trim cigar aligment from the first (match)
function maxTrim1(cigar)
{
	return cigar+0
}
#and remove overlap from first match
function trim1(cigar,ov){
	return (cigar - ov) substr(cigar, index(cigar, "M"))
}

#how much you can trim cigar aligment from the last (match)
function maxTrim2(cigar ,i)
{
	i = length(cigar) - 1
	while (substr(cigar, i, 1) ~ /[0-9]/)
		--i
	
	return substr(cigar, i + 1) + 0
}
#function trim1(cigar, ov    ,s1,s2,le,op){
#	s1 = ov
#	s2 = ov
#	while (1) {
#		le = (cigar + 0)
#		op = substr(cigar, length(le) + 1, 1)
#		if (op == "M") {
#			if (s1 - le <= 0 && s2 - le <= 0) {
#				if (s1 <= 0 && s2 <= 0) 
#					return cigar
#				return (cigar - max(s1, s2)) substr(cigar, length(le) + 1)
#			}
#			s1 -= le
#			s2 -= le
#		} else if (op == "I") {
#			s1 -= le
#		}
#		else if (op == "D") {
#			s2 -= le
#		}
#		else {
#			print "error parsing cigar string"
#			exit
#		}
#		cigar = substr(cigar, length(le) + 2)
#	}
#}

