#samtools mpileup -q 10 -Q 10 -s `cat sorted_bams`|awk -f pileupParser2.awk|awk -f pileup2posterior.awk
#needs mapping.txt
#Part of Lep-MAP3
BEGIN{
	if (limit7 == "") # maximum quality of a read base
		limit7 = 0.001

	map["a"]=1
	map["A"]=1

	map["c"]=2
	map["C"]=2

	map["g"]=3
	map["G"]=3

	map["t"]=4
	map["T"]=4

	"cat mapping.txt"|getline
	if (NF == 0) {
		print "Error: file mapping.txt not found!"  > "/dev/stderr"
		exit
	}
	s = "CHR\tPOS"

	numIndividuals = 0
	for (i = 1; i <= NF; ++i) {
		mapping[i] = $i
		if (!($i in imapping))
			listOrder[++numIndividuals] = $i
		++imapping[$i]
	}
	print "Number of bams = "  NF > "/dev/stderr"
	print "Number of individuals = "  numIndividuals > "/dev/stderr"
	close("cat mapping.txt")


	for (mi = 1; mi <= numIndividuals; ++mi) {
		s = s "\t" listOrder[mi]
	}
	print s
	FS="\t"

	#ascii characters mapping to their quality
	for (i = 33; i <= 127; ++i) {
		if (i == 33) # q == 0
			q = 0.75 # all bases equal probable
		else if (i == 34) # q == 1
			q = 0.5 * (0.75 + 10 ^ (-0.15))  # One base must be more probable than others, e.g. p \in ]10^-0.15, 0.75[
		else
			q = 10 ^ (-0.1 * (i - 33))
	 	quality[sprintf("%c", i)] = q
	}

	for (i = 33; i <= 127; ++i)
		for (j = 33; j <= 127; ++j) {
			qs1 = sprintf("%c", i)
			qs2 = sprintf("%c", j) 
			p = combineQ(qs1, qs2)
			if (p < limit7)
				p = limit7

			logP[qs1 qs2] = log(p / 3)
			logNP[qs1 qs2] = log(1 - p) 
			logNP2[qs1 qs2] = log(p / 3 + (1 - p)) + log(0.5)
		}


	logHalf = log(0.5)
	logTwo = log(2.0)
}

#function phred2p(q    , p)
#{
#	p = ord[q] - 33
#	if (p == 0)
#		return 0.75;
#	if (p == 1)
#		return 0.725;
#	else
#		return 10^(-0.1 * p)
#}


function combineQ(q1, q2     ,ret,p1,p2)
{
	p1 = quality[q1]
	p2 = quality[q2]
	ret = p1 + p2 - p1 * p2
	if (ret > 0.75)
		ret = 0.75
	return ret
}

function myexp(x)
{
	if (x < -100)
		return 0
	else
		return exp(x)
}

{
	delete prob
	for (i = 5; i <= NF; i+=4)
		for (j = 1; j <= length($i); ++j) {
			a = substr($i, j, 1)
			if (a in map) {
				am = map[a]
				individual = mapping[int(i / 4)] 

				#p = combineQ(substr($(i+1), j, 1), substr($(i+2), j, 1))
				#if (p < limit7)
				#	p = limit7
				#logP = log(p / 3)
				#logNP = log(1 - p) 
				#logNP2 = log(p / 3 + (1 - p)) + logHalf
				#faster by using a table...
				qs = substr($(i+1), j, 1) substr($(i+2), j, 1)

				for (k = 1; k <= 4; ++k)
					if (k == am)
						prob[individual, k, k] += logNP[qs]
					else
						prob[individual, k, k] += logP[qs]

				for (k = 1; k < 4; ++k)
					for (l = k + 1; l <= 4; ++l)
						if (k == am || l == am)
							prob[individual, k, l] += logNP2[qs]
						else
							prob[individual, k, l] += logP[qs]
			}
		}
	s = $1 "\t" $2
	for (mi = 1; mi <= numIndividuals; ++mi) {
		m = listOrder[mi]

		maxp = prob[m,1,1] + 0
		for (k = 1; k <= 4; ++k)
			for (l = k; l <= 4; ++l)
				if (prob[m, k, l] + 0 > maxp)
					maxp = prob[m, k, l] + 0

		for (k = 1; k <= 4; ++k)
			for (l = k; l <= 4; ++l) {
				if (k == 1 && l == 1)
					s = s "\t"
				else
					s = s " "
				s = s myexp(prob[m, k, l] - maxp)
			}
	}
	print s
}
