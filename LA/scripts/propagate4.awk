#propagate script for faster version of Lep-Anchor pipeline
#awk -f propagate4.awk pass=1 chr*[0-9].la pass=2 chr*[0-9].la.err|awk -f pickbed.awk - map_extra.bed >map_extra_prop.bed
BEGIN{
        if (scoreDiff == "")
                scoreDiff = 2   #how much score is needed to assign contigs to chrs (2 pacbio reads or 2kb of alignment)

	T="\t"
}
(pass==1 && /^[^#]/){
	contig = $1 T ($7+0) T ($8+0)
	s1 = $11 + $15 / 2
	score1[$5 T contig] =  s1
	s2 = $14 + $15 / 2
	score2[$5 T contig] =  s2
	if (s1 + s2 > 0) {#do not put contigs with zero score to chromosomes...
		if ((contig in chr)) {
			chr[contig] = 0
			print "contig in multiple chromosomes: " contig >"/dev/stderr"
		} else
			chr[contig] = $5
	}
	contigs[contig]
	chrs[$5]
}
#store link scores between contigs...
(pass==2 && NF==6 && $6+0>=1){
	i1=match($3,/[+-]/)
	if (i1>0) {
		c1=$1
		o1=substr($3,i1,1)
		s1=$2
		e1=substr($3,1,i1-1)

		c2=substr($3,i1+1)
		o2=substr($5,length($5))
		s2=$4
		e2=$5+0
		#print c1,s1,e1,o1,c2,s2,e2,o2 T $6
                contig1=c1 T s1 T e1
               	contig2=c2 T s2 T e2
		#same links can be in many files, use FILENAME as hash as well...
		links2[FILENAME T contig1 T o1 T contig2 T o2]+=$6
		if (pair[contig1 T contig2] == "")
			pair[contig1, ++np[contig1]] = contig2
		pair[contig1 T contig2] = 1
		contigs[contig1]
		contigs[contig2]
	}
}
#do the propagation
#(pass==3){
END{
	#get maximum over all filenames...
	for (i in links2) {
		i2 = substr(i, index(i, T) + 1) #strip filename off
		links[i2] = max(links[i2], links2[i]) # and take max...
	}
	delete links2

	#contigs with unique chromosome
	for (c1 in contigs)
		if ((c1 in chr) && chr[c1] > 0)
			output[chr[c1] T c1]

	#first try to put contigs in multiple chromosomes into one chromosome
	for (c1 in contigs)
		if ((c1 in chr) && chr[c1] == 0) {
			maxS = 0
			maxS2 = 0
			maxchr = 0
			for (cc1 in chrs) {
				s = score1[cc1 T c1] + score2[cc1 T c1]
				if (s > maxS) {
					maxS2 = maxS
					maxS = s
					maxchr = cc1
				} else if (s > maxS2)
					maxS2 = s
			
			}
			if (maxS > 0 && maxS - maxS2 >= scoreDiff) {
				chr[c1] = maxchr
				output[maxchr T c1]
				++assigned
			}
		}

	print ("iteration 0: " (assigned + 0) " contigs") >"/dev/stderr"
	assigned = 0

	ors["+"]
	ors["-"]
	#calculate the score of putting each unassigned contig in each chromosome... score1 and score2 are for both ends of a contig
	#could add restriction that score1 and score2 must be computed against different contigs...
	for (c1 in contigs)
		if ((c1 in chr) && chr[c1] > 0) {
			cc1 =  chr[c1]
			n = np[c1]
			for (i = 1; i <= n;++i) {
				c2 = pair[c1, i]
				if (!(c2 in chr) || chr[c2] == 0)
					for (o1 in ors) {
						score1[cc1 T c2] = max(score1[cc1 T c2], links[c1 T o1 T c2 T "+"])
						score2[cc1 T c2] = max(score2[cc1 T c2], links[c1 T o1 T c2 T "-"])
						score1[cc1 T c2] = max(score1[cc1 T c2], links[c2 T "-" T c1 T o1])
						score2[cc1 T c2] = max(score2[cc1 T c2], links[c2 T "+" T c1 T o1])
					}
			}
		}			
	#arrange contigs by score = score1+score2 - second best ...
	maxScore = 0
	for (c1 in contigs)
		if (!(c1 in chr) || chr[c1] == 0) {
			maxS = 0
			maxS2 = 0
			for (cc1 in chrs) {
				s = score1[cc1 T c1] + score2[cc1 T c1]
				if (s > maxS) {
					maxS2 = maxS
					maxS = s
				} else if (s > maxS2)
					maxS2 = s
			}
			s = maxS - maxS2
			if (s >= scoreDiff) {
				scoreContigs[s, ++count[s]] = c1
				if (s > maxScore)
					maxScore = s
			}
		}
	while (maxScore >= scoreDiff) {
		for (;count[maxScore] >= 1; --count[maxScore]) {
			rp = int(1 + rand() * count[maxScore]) #pick random contig with equal score
			c1 = scoreContigs[maxScore, rp]
			scoreContigs[maxScore, rp] = scoreContigs[maxScore, count[maxScore]] 

			if (!(c1 in chr) || chr[c1] == 0) {
				maxS = 0
				maxS2 = 0
				maxchr = 0
				for (cc1 in chrs) {
					s = score1[cc1 T c1] + score2[cc1 T c1]
					if (s > maxS) {
						maxS2 = maxS
						maxS = s
						maxchr = cc1
					} else if (s > maxS2)
						maxS2 = s
			
				}
				if (maxS - maxS2 == maxScore) {# && maxS - maxS2 >= scoreDiff) {
					output[maxchr T c1]
					chr[c1] = maxchr
					++assigned

					#update score1 and score2 for c1
					n = np[c1]
					for (i = 1; i <= n;++i) {
						c2 = pair[c1, i]
						if (!(c2 in chr) || chr[c2] == 0)
							for (o1 in ors) {
								score1[maxchr T c2] = max(score1[maxchr T c2], links[c1 T o1 T c2 T "+"])
								score2[maxchr T c2] = max(score2[maxchr T c2], links[c1 T o1 T c2 T "-"])
								score1[maxchr T c2] = max(score1[maxchr T c2], links[c2 T "-" T c1 T o1])
								score2[maxchr T c2] = max(score2[maxchr T c2], links[c2 T "+" T c1 T o1])
							}
					}
					#update maxScore, scoreContigs[maxScore, count[maxScore]] and count[maxScore]
					for (i = 1; i <= n;++i) {
						c2 = pair[c1, i]
						if (!(c2 in chr) || chr[c2] == 0) {
							maxS = 0
							maxS2 = 0
							for (cc2 in chrs) {
								s = score1[cc2 T c2] + score2[cc2 T c2]
								if (s > maxS) {
									maxS2 = maxS
									maxS = s
	
								} else if (s > maxS2)
									maxS2 = s
			
							}
							s = maxS - maxS2
							if (s >= scoreDiff) {
								if (s < maxScore)
									scoreContigs[s, ++count[s]] = c2
								else {
									--count[maxScore]
									maxScore = s
									scoreContigs[s, ++count[s]] = c2
									++count[s]
								}
							}
						}
					}

				}
			}#if
		}#for

		#find largest maxScore
		while (maxScore >= scoreDiff && count[maxScore]==0)
			--maxScore
#			print maxScore >"/dev/stderr"

	}#while

					
	print ("iteration " ++numIterations ": " assigned " contigs") >"/dev/stderr"
	assigned = 0

	#put the remaining contigs to all chromosomes with at least scoreDiff score...
	for (c1 in contigs)
		if (!(c1 in chr) || chr[c1] == 0)
			for (cc1 in chrs) {
				s = score1[cc1 T c1] + score2[cc1 T c1]
				if (s >= scoreDiff) {
					output[cc1 T c1]
					++assigned
				}
		
			}
	print ("iteration " ++numIterations ": " assigned " contigs (to possibly multiple chromosomes)") >"/dev/stderr"

	do {
		delete output2
		assigned = 0
		for (c1 in contigs)
			if (!(c1 in chr) || chr[c1] == 0)
				for (cc1 in chrs)
					if ((cc1 T c1) in output) { #not in chr but in output
				                n = np[c1]
				                for (i = 1; i <= n;++i) {
				                        c2 = pair[c1, i]
							isOutput = 0
							for (cc2 in chrs)
								if ((cc2 T c2) in output)
									isOutput = 1
				                        if (isOutput==0) { #c2 is not in output
								s1 = 0
								s2 = 0
				                                for (o1 in ors) {
				                                        s1 = max(s1, links[c1 T o1 T c2 T "+"])
				                                        s2 = max(s2, links[c1 T o1 T c2 T "-"])
				                                        s1 = max(s1, links[c2 T "-" T c1 T o1])
				                                        s2 = max(s2, links[c2 T "+" T c1 T o1])
				                                }
								if (s1 + s2 >= scoreDiff) { 
									output2[cc1 T c2] #put c2 in same chrs as c1...
									++assigned
								}
							}
				                }

					}
		for (i in output2) #output2 =>output
			output[i]
		print ("iteration " ++numIterations ": " assigned " contigs (to possibly multiple chromosomes)") >"/dev/stderr"
	} while (assigned > 0)



	for (i in output) {
		$0 = i
		print $2 T $3 T $4 "\t?\t" $1
	}

}#END
function max(a,b)
{
	if (a >= b)
		return a
	return b
}
